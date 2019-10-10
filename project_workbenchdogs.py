import argparse
import copy
import os
import math
import numpy as np

import constants as c
from helpers import divide_into_equal_passes, center_equally_spaced_points_in_range
from helpers import divide_with_clearout_depths
from routerbit import RouterBit
import templates
import toolpaths
from shopbotfile import ShopBotFile
import workpiece

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Drill an array of constant sized holes.')

    parser.add_argument('--bit-diameter', default=0.5, type=float,
                        help='Diameter of the endmill bit being used.')
    parser.add_argument('--bit-clearout-depth', default=2.0+1/32.0, type=float,
                        help='Maximum depth before a clearout operation is performed.')

    parser.add_argument('--hole-diameter', default=0.75, type=float,
                        help='Diameter of the holes being drilled.')

    parser.add_argument('--hole-depth', default=2.0, type=float,
                        help='Depth of the holes being drilled.')

    parser.add_argument('--hole-spacing', default=2.5, type=float,
                        help='Distance between the holes being drilled.')

    parser.add_argument('--ramp-speed', default=0.8, type=float,
                        help='Speed during corners in inches/second. Applies to tight radius helix only.')

    parser.add_argument('--move-speed', default=4.0, type=float,
                        help='Movement speed in inches/sec. Careful! A helical plunge will move at ramp-speed, but large radii with move at this speed.')

    args = parser.parse_args()

    endmill_bit = RouterBit(args.bit_diameter,
                            clearout_depth_inches = args.bit_clearout_depth,
                            pass_depth_inches = 0.25)

    hole_path = toolpaths.HolePath(bit = endmill_bit,
                                   diameter = args.hole_diameter,
                                   depth = args.hole_depth,
                                   direction = c.DIRECTION_CONVENTIONAL)

    multiple_holes_path = toolpaths.MultipleHolesPath(hole_path)
    num_long = 17 # 43.5"
    num_short = 15 #37.5"
    # The workbench has long and short segments. There is a different number of dogs along
    # the two lengths of boards. The table is symmetric about the x-axis, so x=0 is in the
    # center, and the offsets are symmetric about that.
    x_offsets_and_num_points = []


    if 0: # First set of holes
        y_start = 1.5 + args.hole_diameter/2.0

        # Two tracks on board 0
        board_0 = 4.75 + 1/16.0
        x_sum = 0.5 * board_0
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.5 * board_0

        # One track on board 1
        board_1 = (6 + 1/16.0)
        x_sum += 0.5*board_1
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.5*board_1

        # One track on board 2
        board_2 = (6 + 7/8.0)
        x_sum += 0.5*board_2
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.5*board_2

        # One track on board 3
        board_3 = (8.5 - 1/16.0)
        x_sum += 0.5 * board_3
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.5 * board_3

        # One track on board 4
        board_4 = (8.5)
        x_sum += 0.5*board_4
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.5*board_4
    if 1: # second set of holes
        y_start = 0.5 + args.hole_diameter/2.0

        # One tracks down the middle of board 0
        board_0 = 4.75 + 1/16.0
        x_sum = 0.0
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += board_0

        # Two tracks on board 1
        board_1 = (6 + 1/16.0)
        x_sum += 0.25*board_1
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.5*board_1
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.25*board_1

        # Two tracks on board 2
        board_2 = (6 + 7/8.0)
        x_sum += 0.25*board_2
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.5*board_2
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.25*board_2

        # Two tracks on board 3
        board_3 = (8.5 - 1/16.0)
        x_sum += 0.25 * board_3
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.5 * board_3
        x_offsets_and_num_points.append( ( x_sum, y_start, num_short ) )
        x_sum += 0.25 * board_3

        # Two tracks on board 4
        board_4 = (8.5)
        x_sum += 0.25*board_4
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.5*board_4
        x_offsets_and_num_points.append( ( x_sum, y_start, num_long ) )
        x_sum += 0.25*board_4

    print(x_sum)

    mirror_x_offsets = [(-x, y, num) for (x, y, num) in x_offsets_and_num_points[::-1] if x != 0.0]
    x_idx = 0
    num_holes = 0
    for offset_x, y_start, num_y in (mirror_x_offsets + x_offsets_and_num_points):
        direction = 1 if x_idx % 2 == 0 else -1
        y_offsets = [y_start + i*args.hole_spacing for i in range(num_y)]
        x_idx += 1
        num_holes += num_y
        for offset_y in y_offsets[::direction]:
            multiple_holes_path.add_location(offset_x = offset_x, offset_y = offset_y)

    shop_bot_file = ShopBotFile('workbench_holes.sbp')
    shop_bot_file.set_ramps(small_circle_diameter = 0.2,
                            z_move_ramp_speed = args.ramp_speed,
                            z_jog_ramp_speed = args.ramp_speed,
                            xy_move_ramp_speed = args.ramp_speed,
                            xy_jog_ramp_speed = args.ramp_speed)
    shop_bot_file.set_speed(z_move_speed = args.move_speed,
                            z_jog_speed = args.move_speed,
                            xy_move_speed = args.move_speed,
                            xy_jog_speed = args.move_speed)

    shop_bot_file.add_points(multiple_holes_path.get_points())
    shop_bot_file.close()

    points = np.array(multiple_holes_path.get_points())

    width_x = 69.5
    width_y = 43.5
    origin_x = -0.5 * width_x
    origin_y = 0.0
    workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.hole_depth],
                                                   origin=[origin_x, origin_y, 0])

    time_per_hole = 20.0 # seconds (measured)
    path_estimate = workpiece_preview.estimate_time(points, args.ramp_speed)
    empirical_estimate = time_per_hole * num_holes
    print('Emperically estimating {:0.2f} minutes'.format(empirical_estimate / 60.0))
    print('Estimating {:0.2f} minutes based on tool path'.format(path_estimate / 60.0))


    workpiece_preview.plot_birds_eye(points)