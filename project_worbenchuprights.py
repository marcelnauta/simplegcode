import argparse
import copy
import os
import math
import matplotlib.pyplot as plt
import numpy as np

import constants as c
from helpers import divide_into_equal_passes, center_equally_spaced_points_in_range
from helpers import divide_with_clearout_depths
from routerbit import RouterBit
import templates
import toolpaths
from shopbotfile import ShopBotFile
import workpiece

# Effectively a class derived from np.array
def Point2D(x, y):
    return np.array([x, y])

# Intersection of two lines in parametric form
def intersect_lines(point_1, direction_1, point_2, direction_2):
    M = np.array([[-direction_1[0], direction_2[0]],
                  [-direction_1[1], direction_2[1]]])
    X = point_1 - point_2
    parameters = np.dot(np.linalg.inv(M), X)
    return point_1 + parameters[0] * direction_1



def draw_closed_polygon(ax, points, rotation = 0, pivot = [0,0], **kwargs):
    # Subtract pivot
    points[0,:] -= pivot[0]
    points[1,:] -= pivot[1]
    b = np.array(points, copy=True)
    angle = rotation * np.pi / 180.0
    b[0,:] = + np.cos(angle) * points[0,:] + np.sin(angle) * points[1,:]
    b[1,:] = - np.sin(angle) * points[0,:] + np.cos(angle) * points[1,:]
    # Add back pivot
    b[0,:] += pivot[0]
    b[1,:] += pivot[1]
    ax.plot(b[0,:], b[1,:], **kwargs)

def draw_box(ax, origin, size, **kwargs):
    a = np.array([[origin[0], origin[0]+size[0], origin[0]+size[0], origin[0]        , origin[0] ],
                  [origin[1], origin[1]        , origin[1]+size[1], origin[1]+size[1], origin[1] ]])
    draw_closed_polygon(ax, a, **kwargs)

def draw_circle(ax, center, diameter, **kwargs):
    N = 64
    a = np.zeros((2,N))
    theta = np.linspace(0, 2*np.pi, N)
    a[0,:] = center[0] + 0.5 * diameter * np.cos(theta)
    a[1,:] = center[1] + 0.5 * diameter * np.sin(theta)
    draw_closed_polygon(ax, a, **kwargs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Drill an array of constant sized holes.')

    parser.add_argument('--bit-diameter', default=0.5, type=float,
                        help='Diameter of the endmill bit being used.')
    parser.add_argument('--bit-clearout-depth', default=2.0+1/32.0, type=float,
                        help='Maximum depth before a clearout operation is performed.')

    parser.add_argument('--axle-hole-diameter', default=1.0, type=float,
                        help='Diameter of the hinge axle.')

    parser.add_argument('--thickness', default=1.0 + 7/8.0, type=float,
                        help='Thickness of the workpiece')

    parser.add_argument('--cart-width', default=20.0, type=float,
                        help='Width of the cart')
    parser.add_argument('--cart-thickness', default=1.0 + 7/8.0, type=float,
                        help='Thickness of the cart')
    parser.add_argument('--cart-height-at-bottom', default=7.0, type=float,
                        help='Distance between the floor and the bottom of the cart')

    parser.add_argument('--pivot-fraction-of-cart', default=0.45, type=float,
                        help='Set the pivot location in x as a fraction of the cart width')
    parser.add_argument('--pivot-fraction-of-upright', default=0.85, type=float,
                        help='Set the pivot location in x as a fraction of the width of the upright')
    parser.add_argument('--pivot-fraction-of-table', default=0.55, type=float,
                        help='Set the pivot location in x as a fraction of the width of the table')
    parser.add_argument('--pivot-to-table-distance', default=3.5, type=float,
                        help='Distance between the center of the pivot and the bottom of the table (min of sqrt(2)/2 * flange_width)')
    parser.add_argument('--pivot-max-angle', default=78.0, type=float,
                        help='Maximum angle of the pivot in degrees')

    parser.add_argument('--bearing-flange-size', default=3 + 5/8.0, type=float,
                        help='Outside to outside size of the flange bearing')

    parser.add_argument('--table-height', default=40.0, type=float,
                        help='Desired height of the top of the table')
    parser.add_argument('--table-thickness', default=1 + 7/8.0, type=float,
                        help='Thickness of the table')
    parser.add_argument('--table-width', default=42.0, type=float,
                        help='Wdith of the table')

    parser.add_argument('--vertical-board-width', default=4.0 + 3/16.0, type=float,
                        help='Width of the straight vertical board')
    parser.add_argument('--angled-board-width', default=7.0, type=float,
                        help='Width of the angled vertical board')
    parser.add_argument('--angled-board-length', default=38.0, type=float,
                        help='Uncut length of the angled vertical board')
    parser.add_argument('--brace-board-width', default=5.0 + 5/8.0, type=float,
                        help='Width of the thinner horizontal board')
    parser.add_argument('--brace-board-thickness', default=7.0/8.0 + 1.0/32.0, type=float,
                        help='Thickness of the thinner horizontal board')

    parser.add_argument('--lap-ramp-speed', default=0.4, type=float,
                        help='Movement speed in inches/sec along tight curves within lap joint')

    parser.add_argument('--lap-move-speed', default=3.5, type=float,
                        help='Movement speed in inches/sec along straight paths of lap joint.')

    parser.add_argument('--ramp-speed', default=0.8, type=float,
                        help='Speed during corners in inches/second. Applies to tight radius helix only.')

    parser.add_argument('--move-speed', default=1.0, type=float,
                        help='Movement speed in inches/sec. Careful! A helical plunge will move at ramp-speed, but large radii with move at this speed.')


    args = parser.parse_args()

    # Draw a diagram of the project
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1,1,1, aspect='equal')

    # Compute some derived constants
    upright_length = args.table_height - args.table_thickness - args.cart_height_at_bottom
    pivot_x = (-0.5 + args.pivot_fraction_of_cart) * args.cart_width
    upright_start_x = pivot_x - args.pivot_fraction_of_upright * args.vertical_board_width
    pivot_y = args.table_height - args.table_thickness - args.pivot_to_table_distance

    # The corners of the two uprights are found using the following method:
    #    a) Find a point by starting at the pivot point, and going out by +/- thickness/2 at the angle perpindicular to the board.
    #    b) Find the equation for the line going through that point, parallel to the board.
    #    c) Find the point that intersects the bottom of the table, when lying flat.
    #    d) Find the point that intersects the bottom of the cart, which always lies flat.
    for board_angle_degrees in [args.pivot_max_angle]:
        upright_points = []
        for direction in [+1, -1]:
            half_width = direction * 0.5 * args.angled_board_width
            board_angle = board_angle_degrees * np.pi / 180.0
            point_on_line = Point2D(pivot_x, pivot_y) + half_width * Point2D(math.sin(board_angle), math.cos(board_angle))
            board_direction = Point2D(math.cos(board_angle), -math.sin(board_angle))
            point_at_table = intersect_lines(point_on_line, board_direction,
                                            Point2D(0,args.table_height - args.table_thickness),
                                            Point2D(1,0))
            point_at_cart = intersect_lines(point_on_line, board_direction,
                                            Point2D(0,args.cart_height_at_bottom),
                                            Point2D(1,0))
            if direction > 0:
                upright_points.extend([point_at_table, point_at_cart])
            else:
                upright_points.extend([point_at_cart, point_at_table])
        draw_closed_polygon(ax, np.transpose(upright_points), color = 'k')

    # Cart base
    draw_box(ax, [-args.cart_width/2.0, args.cart_height_at_bottom], [args.cart_width, args.cart_thickness])
    # Upright
    draw_box(ax, [-args.cart_width/2.0, args.cart_height_at_bottom], [args.vertical_board_width, upright_length])
    # Brace
    draw_box(ax, [-args.cart_width/2.0, args.table_height - args.table_thickness], [args.cart_width, -args.brace_board_width])
    # Table top laying at several angles
    for angle in [0, args.pivot_max_angle / 3.0, 2.0 * args.pivot_max_angle / 3.0, args.pivot_max_angle]:
        draw_box(ax, [-args.pivot_fraction_of_table*args.table_width, args.table_height], [args.table_width, -args.table_thickness],
                 rotation = angle, pivot = [pivot_x, pivot_y], color = 'g')
    # "Wheels"
    draw_box(ax, [-args.cart_width/2.0, 0], [args.cart_width, args.cart_height_at_bottom])
    # Flange bearing mounted at 45 degrees
    draw_box(ax, [pivot_x - args.bearing_flange_size/2.0, pivot_y - args.bearing_flange_size/2.0], [args.bearing_flange_size, args.bearing_flange_size],
             rotation = args.pivot_max_angle, pivot = [pivot_x, pivot_y])
    draw_circle(ax, [pivot_x, pivot_y], args.axle_hole_diameter, color = 'k')

    ## SETUP CUTTING PROFILES

    endmill_bit = RouterBit(args.bit_diameter, clearout_depth_inches = args.bit_clearout_depth, pass_depth_inches = 0.25)
    # Step 1) Cut out an angled lap joint between the angled board and the horizontal brace
    # Step 2) Manually glue the lap joint together
    # Step 3) Cut out the hole for the axle and trim the angles of the board
    for upright_name in ['left', 'right']:
        angle_length_multiplier = 1.0 / math.sin(args.pivot_max_angle * np.pi / 180.0)
        dato_multiplier = (angle_length_multiplier - 1) * 2.0 + 1.0 + 0.2
        board_cutoff = 1.0
        rotation = (args.pivot_max_angle - 90) if upright_name == 'left' else (90 - args.pivot_max_angle)

        upper_lap_path = toolpaths.DatoPath(bit = endmill_bit,
                                            dato_width = args.brace_board_width,
                                            board_width = args.angled_board_width * dato_multiplier,
                                            depth = args.brace_board_thickness,
                                            direction = c.DIRECTION_CONVENTIONAL)
        upper_lap_path.add_transform(toolpaths.Rotate2D(rotation))
        upper_lap_height_wrt_cart = args.table_height - args.table_thickness - args.brace_board_width / 2.0 - args.cart_height_at_bottom + board_cutoff
        upper_lap_path.add_transform(toolpaths.Shift2D( args.angled_board_width/2.0 , upper_lap_height_wrt_cart * angle_length_multiplier))

        lower_lap_path = toolpaths.DatoPath(bit = endmill_bit,
                                            dato_width = args.cart_thickness,
                                            board_width = args.angled_board_width * dato_multiplier,
                                            depth = args.brace_board_thickness,
                                            direction = c.DIRECTION_CONVENTIONAL)
        lower_lap_path.add_transform(toolpaths.Rotate2D(rotation))
        lower_lap_height_wrt_cart = args.cart_thickness/2.0 + board_cutoff
        lower_lap_path.add_transform(toolpaths.Shift2D( args.angled_board_width/2.0 , lower_lap_height_wrt_cart * angle_length_multiplier))

        axle_hole_path = toolpaths.HolePath(bit = endmill_bit,
                                            diameter = args.axle_hole_diameter,
                                            depth = args.thickness,
                                            direction = c.DIRECTION_CONVENTIONAL)
        axle_hole_height_wrt_cart = pivot_y - args.cart_height_at_bottom + board_cutoff
        axle_hole_path.add_transform(toolpaths.Shift2D(args.angled_board_width/2.0, axle_hole_height_wrt_cart * angle_length_multiplier))

        lap_shop_bot_file = ShopBotFile(upright_name + '_angled_upright_lap_joints.sbp')
        lap_shop_bot_file.set_ramps(small_circle_diameter = 0.2,
                                    z_move_ramp_speed = args.lap_ramp_speed,
                                    z_jog_ramp_speed = args.lap_ramp_speed,
                                    xy_move_ramp_speed = args.lap_ramp_speed,
                                    xy_jog_ramp_speed = args.lap_ramp_speed)
        lap_shop_bot_file.set_speed(z_move_speed = args.lap_move_speed,
                                    z_jog_speed = args.lap_move_speed,
                                    xy_move_speed = args.lap_move_speed,
                                    xy_jog_speed = args.lap_move_speed)

        hole_shop_bot_file = ShopBotFile(upright_name + '_angled_upright_hole.sbp')
        hole_shop_bot_file.set_ramps(small_circle_diameter = 0.2,
                                     z_move_ramp_speed = args.ramp_speed,
                                     z_jog_ramp_speed = args.ramp_speed,
                                     xy_move_ramp_speed = args.ramp_speed,
                                     xy_jog_ramp_speed = args.ramp_speed)
        hole_shop_bot_file.set_speed(z_move_speed = args.move_speed,
                                     z_jog_speed = args.move_speed,
                                     xy_move_speed = args.move_speed,
                                    xy_jog_speed = args.move_speed)

        lap_shop_bot_file.add_points(lower_lap_path.get_points())
        lap_shop_bot_file.add_points(upper_lap_path.get_points())
        hole_shop_bot_file.add_points(axle_hole_path.get_points())
        lap_shop_bot_file.close()
        hole_shop_bot_file.close()

        upper_lap_points = np.array(upper_lap_path.get_points())
        lower_lap_points = np.array(lower_lap_path.get_points())
        axle_hole_points = np.array(axle_hole_path.get_points())

        width_x = args.angled_board_width
        origin_x = 0.0
        width_y = args.angled_board_length
        origin_y = 0.0
        workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.thickness],
                                                    origin=[origin_x, origin_y, 0])
        workpiece_preview.plot_birds_eye(upper_lap_points, extras = [lower_lap_points, axle_hole_points])
        #workpiece_preview.plot_three_axis(axle_hole_points, draw_bounding_box = False)

    # Step 4) Make a bracket to attach the flange bearing to the table top.
    #    14" wide
    #    3 + 11/16" flange
    #    taper between two thickness with a cos**2 taper
    #    save material by cutting the boards out with a bandsaw first

    '''
    hole_path = toolpaths.HolePath(bit = endmill_bit,
                                   diameter = args.hole_diameter,
                                   depth = args.hole_depth,
                                   direction = c.DIRECTION_CONVENTIONAL)

    multiple_holes_path = toolpaths.MultipleHolesPath(hole_path)
    x_offsets = center_equally_spaced_points_in_range(0, args.x_length, args.x_spacing, args.x_min_edge)
    y_offsets = center_equally_spaced_points_in_range(0, args.y_length, args.y_spacing, args.y_min_edge)
    print('X locations: {}'.format(x_offsets))
    print('Y locations: {}'.format(y_offsets))
    for x_idx, offset_x in enumerate(x_offsets):
        direction = 1 if x_idx % 2 == 0 else -1
        for offset_y in y_offsets[::direction]:
            multiple_holes_path.add_location(offset_x = offset_x, offset_y = offset_y)

    shop_bot_file = ShopBotFile('custom_size_hole.sbp')
    shop_bot_file.set_ramps(small_circle_diameter = 0.2, xy_move_ramp_speed = 0.8)

    shop_bot_file.add_points(multiple_holes_path.get_points())
    shop_bot_file.close()

    points = np.array(multiple_holes_path.get_points())

    width_x = args.x_length
    origin_x = 0.0
    if width_x == 0.0:
        width_x = 2.0 * args.hole_diameter
        origin_x = -0.5 * width_x

    width_y = args.y_length
    origin_y = 0.0
    if width_y == 0.0:
        width_y = 2.0 * args.hole_diameter
        origin_y = -0.5 * width_y
    workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.hole_depth],
                                                   origin=[origin_x, origin_y, 0])
    workpiece_preview.plot_three_axis(points)
    '''