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

from variablelocations import finger_factory

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Drill an array of constant sized holes.')

    parser.add_argument('--bit-diameter', default=0.25, type=float,
                        help='Diameter of the endmill bit being used.')
    parser.add_argument('--bit-clearout-depth', default=0.5+1/32.0, type=float,
                        help='Maximum depth before a clearout operation is performed.')

    parser.add_argument('--total-length', default=24.0, type=float,
                        help='Thickness of each of the slats')

    parser.add_argument('--slat-thickness', default=0.375, type=float,
                        help='Thickness of each of the slats')

    parser.add_argument('--min-slat-width', default=5/8.0, type=float,
                        help='Min width of the slats')

    parser.add_argument('--max-slat-width', default=12/8.0, type=float,
                        help='Max width of the slats')

    parser.add_argument('--hole-depth', default=1.25, type=float,
                        help='Depth of the holes being drilled.')

    parser.add_argument('--slat-spacing', default=2.25, type=float,
                        help='Distance between the holes being drilled.')

    parser.add_argument('--ramp-speed', default=0.8, type=float,
                        help='Speed during corners in inches/second. Applies to tight radius helix only.')

    parser.add_argument('--move-speed', default=1.25, type=float,
                        help='Movement speed in inches/sec. Careful! A helical plunge will move at ramp-speed, but large radii with move at this speed.')

    args = parser.parse_args()

    endmill_bit = RouterBit(args.bit_diameter,
                            clearout_depth_inches = args.bit_clearout_depth,
                            pass_depth_inches = 5/32.0)

    entry = finger_factory('SquishMulti', 'finger_style', points_from_json=True)
    if 1:
        entry_json = '''
{
    "finger_style": {
        "min_size": 2.5,
        "min_inset": 8.0,
        "max_inset": 8.0,
        "max_size": 6.5,
        "periods": 2.0
    }
}
'''
    entry_args = {'finger_style': {
            'min_size': args.slat_spacing + args.min_slat_width,
            'max_size': args.slat_spacing + args.max_slat_width,
            'min_inset': 0.0,
            'max_inset': 0.0,
            'periods': 1.0,
            }}
    entry.from_property_tree(entry_args)
    entry.length = args.total_length
    locations = [0.0, ] + entry.get_locations() + [args.total_length, ]
    print(locations)

    points_list = []
    for shifted_i, y_loc in enumerate(locations[1:-1]):
        i = shifted_i + 1
        width_y = 0.5*(locations[i+1] - locations[i-1] - 2*args.slat_spacing)
        print(width_y)
        hole_path = toolpaths.SquareHolePath(bit = endmill_bit,
                                             width_x = args.slat_thickness,
                                             width_y = width_y,
                                             depth = args.hole_depth,
                                             direction = c.DIRECTION_CONVENTIONAL,
                                             offset_y = y_loc,
                                             rotation = -90.0)
        points_list.extend(hole_path.get_points())


    shop_bot_file = ShopBotFile('crib_slats_{}.sbp'.format(int(args.total_length)))
    shop_bot_file.set_ramps(small_circle_diameter = 0.2,
                            z_move_ramp_speed = args.ramp_speed,
                            z_jog_ramp_speed = args.ramp_speed,
                            xy_move_ramp_speed = args.ramp_speed,
                            xy_jog_ramp_speed = args.ramp_speed)
    shop_bot_file.set_speed(z_move_speed = args.move_speed,
                            z_jog_speed = args.move_speed,
                            xy_move_speed = args.move_speed,
                            xy_jog_speed = args.move_speed)

    shop_bot_file.add_points(points_list)
    shop_bot_file.close()

    points = np.array(points_list)

    width_x = args.total_length
    width_y = 7/8.0
    origin_x = 0.0
    origin_y = -0.5 * width_y
    workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.hole_depth],
                                                   origin=[origin_x, origin_y, 0])

    path_estimate = workpiece_preview.estimate_time(points, args.ramp_speed)
    print('Estimating {:0.2f} minutes based on tool path'.format(path_estimate / 60.0))

    workpiece_preview.plot_birds_eye(points)
    workpiece_preview.plot_three_axis(points)