import copy
import math
import numpy as np

import constants as c
from helpers import divide_into_equal_passes, divide_with_clearout_depths

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

class Rotate2D(PathTransformABC):
    def __init__(self, angle_ccw, degrees = True):
        scalar = math.pi / 180.0 if degrees else 1.0
        self.angle_ccw = angle_ccw * scalar

    def apply_transform(self, point):
        rv = copy.deepcopy(point)
        rv[c.AXIS_X] = math.cos(self.angle_ccw) * point[c.AXIS_X] - math.sin(self.angle_ccw) * point[c.AXIS_Y]
        rv[c.AXIS_Y] = math.sin(self.angle_ccw) * point[c.AXIS_X] + math.cos(self.angle_ccw) * point[c.AXIS_Y]
        return rv

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
        num_rings = max(1, math.ceil((max_radius - min_radius) / self.bit.pass_horiz) + 1)

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
        # Only raise slightly above the bottom so that we can move at jog speed up out of the hole.
        last_known_radius = points[-1][1]
        for theta in full_circle_angles:
            depth_scalar = math.pow(math.cos(0.25 * theta), 0.25)
            withdraw_depth = max(depth - self.bit.pass_depth, 0.0)
            radial_scalar = 1.0 - (theta / (2*math.pi))**2
            this_radius = radial_scalar * last_known_radius
            points.append([this_radius * math.sin(theta), this_radius * math.cos(theta), -(depth_scalar*depth + (1-depth_scalar)*withdraw_depth)])

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


class DatoPath(ToolPathABC):
    def __init__(self, bit, dato_width, board_width, depth,
                 direction = c.DIRECTION_CONVENTIONAL,
                 offset_x = 0.0, offset_y = 0.0):
        ToolPathABC.__init__(self)
        self.bit = bit
        self.diameter = dato_width
        self.dato_width = dato_width
        self.board_width = board_width
        self.depth = depth
        if direction == c.DIRECTION_CLIMB:
            self.add_transform(Mirror(c.AXIS_Y))
        if offset_x != 0.0 or offset_y != 0.0:
            self.add_transform(Shift2D(offset_x, offset_y))

    def get_untransformed_points(self):
        # Untransformed dato assumes the board is lying in the y direction with the cut in the x direction
        # The origin is the center of the dato in the center of the board.
        min_x = -self.board_width / 2.0 - self.bit.radius # outside workpiece
        max_x = -min_x
        min_y = -self.dato_width / 2.0 + self.bit.radius # inside slot size
        max_y = -min_y
        points = [[min_x, min_y, c.SAFE_HEIGHT]]

        depths = divide_with_clearout_depths(self.depth, self.bit.pass_depth, self.bit.clearout_depth)
        for depth in depths:
            x_off = 0.0
            for y in divide_into_equal_passes(max_y, self.bit.pass_horiz/2.0, self.bit.pass_horiz)[::-1]:
                points.append([min_x - x_off, -y, -depth])
                points.append([min_x - x_off,  y, -depth])
                points.append([max_x + x_off,  y, -depth])
                points.append([max_x + x_off, -y, -depth])
                points.append([min_x - x_off, -y, -depth])
                #x_off -= self.bit.pass_horiz

        points.append([points[-1][0] ,points[-1][1], c.SAFE_HEIGHT])
        return points


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