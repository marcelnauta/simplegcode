import copy
import os
import math
import numpy as np

import constants as c
from helpers import divide_into_equal_passes, center_equally_spaced_points_in_range
from helpers import divide_with_clearout_depths
import templates
import workpiece

class ShopBotFile(object):
    def __init__(self, filename = 'tmp.sbp', xy_speed = 1.5, z_speed = 1.5, a_speed = 360.0, b_speed = 360.0):
        self.outfile = open(filename, 'w')
        self._write(templates.get_header())
        self.movement_speeds = c.MovementSpeeds()
        self.ramp_speeds = c.RampSpeeds()
        self.set_speed()
        self.set_ramps(small_circle_diameter = 0.2,
                       xy_move_ramp_speed = 0.8)

    def _write(self, out_str):
        self.outfile.write(out_str)
        #print(out_str)

    def close(self):
        self._write(templates.get_footer())
        self.outfile.close()

    def set_speed(self, **kwargs):
        self.movement_speeds.set(**kwargs)
        self._write(self.movement_speeds.get_command())

    def set_ramps(self, **kwargs):
        self.ramp_speeds.set(**kwargs)
        self._write(self.ramp_speeds.get_command())

    def move_to(self, point, command_prefix = 'M'):
        if len(point) == 4: # At Maker Labs, the CNC is 4 axis, but the 4th axis acts as the 5th
            point = (point[0], point[1], point[2], 0.0, point[4])
        point_as_str = ['{0:.6f}'.format(x) for x in point]
        command_tag = command_prefix + '{0:d}'.format(len(point))
        self._write(','.join([command_tag] + point_as_str) + '\n')

    def jog_to(self, point):
        self.move_to(point, command_prefix = 'J')

    def add_points(self, point_list):
        [self.move_to(point) for point in point_list]

    def pause(self, time):
        self._write('PAUSE {}\n'.format(time))

class RouterBit(object):
    def __init__(self, diameter):
        self.diameter = diameter
        self.radius = diameter / 2.0
        self.pass_depth = 3/16.0
        self.clearout_depth = 2.01
        self.pass_horiz = 0.3 * diameter

    def get_radius(height_wrt_bottom):
        return self.radius

class PathTransformABC(object):
    def apply_transform(point):
        raise NotImplementedError('apply_transform should be pure virtual')

    def apply_transform_to_list(self, points):
        return [self.apply_transform(copy.deepcopy(point)) for point in points]

class Shift2D(PathTransformABC):
    def __init__(self, offset_x, offset_y):
        self.offset_x = offset_x
        self.offset_y = offset_y

    def apply_transform(self, point):
        point[c.AXIS_X] = point[c.AXIS_X] + self.offset_x
        point[c.AXIS_Y] = point[c.AXIS_Y] + self.offset_y
        return point

class Mirror(PathTransformABC):
    def __init__(self, axis):
        self.axis = axis

    def apply_transform(self, point):
        point[self.axis] = -point[self.axis]
        return point

class ToolPathABC(object):
    def __init__(self):
        self.paths = []
        self.transforms = []

    def get_untransformed_points(self):
        return []

    def get_points(self):
        points = []
        for path in self.paths:
            points.extend(path.get_points())
        points.extend(self.get_untransformed_points())
        for transform in self.transforms:
            points = transform.apply_transform_to_list(points)
        return points

    def add_transform(self, transform):
        self.transforms.append(transform)

class MoveToPath(ToolPathABC):
    def __init__(self, point):
        ToolPathABC.__init__(self)
        self.point = point

    def get_untransformed_points(self):
        return [self.point]

class HolePath(ToolPathABC):
    def __init__(self, bit, diameter, depth, direction = c.DIRECTION_CONVENTIONAL,
                 offset_x = 0.0, offset_y = 0.0):
        ToolPathABC.__init__(self)
        self.bit = bit
        self.diameter = diameter
        self.depth = depth
        if direction == c.DIRECTION_CLIMB:
            self.add_transform(Mirror(c.AXIS_Y))
        if offset_x != 0.0 or offset_y != 0.0:
            self.add_transform(Shift2D(offset_x, offset_y))

    def get_untransformed_points(self):
        max_radius = max(0.0, self.diameter / 2.0 - self.bit.radius)
        min_radius = min(max_radius, 0.8 * self.bit.radius)
        num_rings = max(1, math.ceil((max_radius - min_radius) / self.bit.pass_horiz))

        full_circle_angles = np.linspace(0, 2*math.pi, 65, endpoint=True)

        points = [[0,min_radius,c.SAFE_HEIGHT]]

        # Cork screw down into the work piece in case the material does not start precisely at
        # z = 0 and to prevent stopping on the material surface.
        descent_heights = divide_into_equal_passes(c.SAFE_HEIGHT, 0.0, self.bit.pass_depth, include_first = True)
        num_descent_steps = len(descent_heights)
        for height_idx in range(1, num_descent_steps):
            for theta in full_circle_angles:
                height_scalar = math.sin(0.25 * theta)**2
                this_height = height_scalar * descent_heights[height_idx] + (1.0-height_scalar) * descent_heights[height_idx-1]
                points.append([min_radius * math.sin(theta), min_radius * math.cos(theta), this_height])

        # Without ever clearing sawdust out, the hole would be cut in the following passes
        depths = divide_with_clearout_depths(self.depth, self.bit.pass_depth, self.bit.clearout_depth)

        if max_radius == 0.0:
            for depth in depths:
                points.append([0,0,0])
                points.append([0,0,-depth])

        for depth_idx, depth in enumerate(depths):
            new_depth = True
            radii = np.linspace(min_radius, max_radius, num_rings, endpoint=True)
            # Ramp down in depth over one ring
            for theta in full_circle_angles:
                depth_scalar = math.sin(0.25 * theta)**2
                depth_start = depths[depth_idx-1] if depth_idx != 0 else 0.0
                this_depth = depth_scalar * depths[depth_idx] + (1.0-depth_scalar) * depth_start
                points.append([min_radius * math.sin(theta), min_radius * math.cos(theta), -this_depth])
            # Do one full ring for the first radius at the current depth
            if min_radius != max_radius or depth_idx == len(depths)-1:
                for theta in full_circle_angles:
                    points.append([min_radius * math.sin(theta), min_radius * math.cos(theta), -depth])

            if min_radius != max_radius:
                # Spiral out to the outer radius
                for radius_idx in range(len(radii)-1):
                    inner_radius = radii[radius_idx]
                    outer_radius = radii[radius_idx+1]
                    for theta in full_circle_angles:
                        # taper the spiral profile so that it departs and arrives tangential to the rings
                        taper_scalar = theta / (2*math.pi)
                        taper_width = 0.2
                        if taper_scalar < taper_width:# and radius_idx == 0:
                            radius_scalar = math.sin(taper_scalar / taper_width * 0.5 * np.pi)**2 * taper_width
                        elif taper_scalar > 1.0 - taper_width:# and radius_idx == len(radii)-2:
                            tmp_scalar = 1.0 - taper_scalar
                            radius_scalar = 1.0 - math.sin(tmp_scalar / taper_width * 0.5 * np.pi)**2 * tmp_scalar
                        else:
                            radius_scalar = taper_scalar
                        this_radius = radius_scalar * outer_radius + (1.0-radius_scalar) * inner_radius
                        points.append([this_radius * math.sin(theta), this_radius * math.cos(theta), -depth])

                # Do one full ring for the last radius at the current depth
                for theta in full_circle_angles:
                    points.append([max_radius * math.sin(theta), max_radius * math.cos(theta), -depth])

                # Spiral back to the inner radius in a single circle at the current depth
                # Note, spiralling back prevents stopping to make a 90 degree turn. The
                # radius of the spiral is tapered to make the changes in direction even smoother
                for theta in full_circle_angles:
                    radius_scalar = math.sin(theta / 4.0)**2
                    this_radius = radius_scalar * min_radius + (1.0-radius_scalar) * max_radius
                    points.append([this_radius * math.sin(theta), this_radius * math.cos(theta), -depth])

            #for radius in np.linspace(min_radius, max_radius, num_rings, endpoint=True):
            #    for theta in full_circle_angles:
            #        points.append([radius * math.sin(theta), radius * math.cos(theta), -depth])

        # Spiral back into the center from the last used depth. Increase height very non-linearly
        last_known_radius = points[-1][1]
        for theta in full_circle_angles:
            depth_scalar = math.pow(math.cos(0.25 * theta), 0.25)
            radial_scalar = 1.0 - (theta / (2*math.pi))**2
            this_radius = radial_scalar * last_known_radius
            points.append([this_radius * math.sin(theta), this_radius * math.cos(theta), -depth_scalar*depth])

        points.append([0,0,c.SAFE_HEIGHT])
        return points

class MultipleHolesPath(ToolPathABC):
    def __init__(self, hole_path, xy_locations = []):
        ToolPathABC.__init__(self)
        self.hole_path = hole_path
        for xy_location in xy_locations:
            self.add_location(xy_location[0], xy_location[1])

    def add_location(self, offset_x = 0.0, offset_y = 0.0):
        new_hole_path = copy.deepcopy(self.hole_path)
        new_hole_path.add_transform(Shift2D(offset_x, offset_y))
        self.paths.append(new_hole_path)


def make_dovetail_domino(width = 5/8.0, half_height = 7/16.0, length = 8.0,
                         pass_width = 1/16.0, bit_width = 1/2.0, num_passes = 2):
    point_list = []
    start_point = [0, 0, c.SAFE_HEIGHT]
    half_length = length/2.0

    point_list.append(start_point)

    pass_offsets = [i * pass_width for i in range(num_passes, -1, -1)]
    for pass_offset in pass_offsets:
        this_pass_width = width/2.0 + bit_width/2.0 + pass_offset
        # Drop the bit next to the work piece
        plunge_width = this_pass_width + bit_width/2.0
        point_list.append([half_length, -plunge_width, c.SAFE_HEIGHT])
        point_list.append([half_length, -plunge_width, -half_height])
        # Approach the work piece
        point_list.append([half_length, -this_pass_width, -half_height])
        # Run down the length of it
        point_list.append([-half_length, -this_pass_width, -half_height])
        # Pull away and raise up
        point_list.append([-half_length, -plunge_width, -half_height])
        point_list.append([-half_length, -plunge_width, c.SAFE_HEIGHT])

        # Drop the bit next to the work piece on the other side
        plunge_width = this_pass_width + bit_width/2.0
        point_list.append([-half_length, plunge_width, c.SAFE_HEIGHT])
        point_list.append([-half_length, plunge_width, -half_height])
        # Approach the work piece
        point_list.append([-half_length, this_pass_width, -half_height])
        # Run up the length of it
        point_list.append([half_length, this_pass_width, -half_height])
        # Pull away and raise up
        point_list.append([half_length, plunge_width, -half_height])
        point_list.append([half_length, plunge_width, c.SAFE_HEIGHT])
    point_list.append(start_point)
    shop_bot_file = ShopBotFile('dovetail_dominoes.sbp')
    shop_bot_file.add_points(point_list)
    shop_bot_file.close()

    # Make a holder to consistently align the dominoes.
    # The bottom will define the new z=0
    point_list = []
    point_list.append(start_point)
    pass_heights = np.linspace(0.25, 0.0, num = 10, endpoint=True)
    widths =  np.arange(width /2.0 - bit_width/2.0, 0.0, -bit_width/2.0)

    first_pass = True

    for pass_height in pass_heights:
        this_length = half_length + bit_width/2.
        for this_width in  widths:
            if first_pass:
                first_pass = False
                point_list.append([ this_length, -this_width, c.SAFE_HEIGHT])
            point_list.append([ this_length, -this_width, pass_height])
            point_list.append([-this_length, -this_width, pass_height])
            point_list.append([-this_length,  this_width, pass_height])
            point_list.append([ this_length,  this_width, pass_height])
            point_list.append([ this_length, -this_width, pass_height])
    point_list.append([ this_length,  this_width, c.SAFE_HEIGHT])

    point_list.append(start_point)
    shop_bot_file = ShopBotFile('dovetail_dominoes_holder.sbp')
    shop_bot_file.add_points(point_list)
    shop_bot_file.close()




if __name__ == '__main__':
    #make_dovetail_domino()
    dovetail_bit = RouterBit(0.5)
    endmill_bit = RouterBit(0.5)

    move_to_path = MoveToPath([12.0, 2.0, 0.0, math.pi])
    shift_2d = Shift2D(1.0, 1.0)
    move_to_path.add_transform(shift_2d)

    hole_path = HolePath(bit = endmill_bit, diameter = 0.75, depth = 1.55,
                         direction = c.DIRECTION_CONVENTIONAL)


    multiple_holes_path = MultipleHolesPath(hole_path)
    x_offsets = center_equally_spaced_points_in_range(0, 42, 2.5, 2.0)
    y_offsets = [0.0] #center_equally_spaced_points_in_range(0, 42, 2.5, 1.0)
    for x_idx, offset_x in enumerate(x_offsets):
        direction = 1 if x_idx % 2 == 0 else -1
        for offset_y in y_offsets[::direction]: #y_offsets[::direction]:
            multiple_holes_path.add_location(offset_x = offset_x, offset_y = offset_y)

    shop_bot_file = ShopBotFile('custom_size_hole.sbp')

    shop_bot_file.add_points(multiple_holes_path.get_points())
    shop_bot_file.close()

    points = np.array(multiple_holes_path.get_points())

    workpiece_preview = workpiece.WorkpiecePreview([2, 2, 2], origin=[-1, -1, 0])
    workpiece_preview.plot_three_axis(points)

    print(points)