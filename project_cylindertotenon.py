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

    parser.add_argument('--tenon-size', default=1.75, type=float,
                        help='Diameter of the holes being drilled.')

    parser.add_argument('--cylinder-size', default=2.25, type=float,
                        help='Depth of the holes being drilled.')

    parser.add_argument('--ramp-speed', default=0.4, type=float,
                        help='Speed during corners in inches/second. Applies to tight radius helix only.')

    parser.add_argument('--move-speed', default=2.0, type=float,
                        help='Movement speed in inches/sec. Careful! A helical plunge will move at ramp-speed, but large radii with move at this speed.')
    
    parser.add_argument('--rot-ramp-speed', default=40.0, type=float,
                        help='')
    parser.add_argument('--rot-move-speed', default=400.0, type=float,
                        help='')

    args = parser.parse_args()

    endmill_bit = RouterBit(args.bit_diameter, clearout_depth_inches = args.bit_clearout_depth)

    #hole_path = toolpaths.HolePath(bit = endmill_bit,
    #                               diameter = args.hole_diameter,
    #                               depth = args.hole_depth,
    #                               direction = c.DIRECTION_CONVENTIONAL)


    shop_bot_file = ShopBotFile('cylinder_to_tenon.sbp')
    shop_bot_file.set_ramps(small_circle_diameter = 0.2,
                            z_move_ramp_speed = args.ramp_speed,
                            z_jog_ramp_speed = args.ramp_speed,
                            xy_move_ramp_speed = args.ramp_speed,
                            xy_jog_ramp_speed = args.ramp_speed,
                            b_move_ramp_speed = args.rot_ramp_speed,
                            b_jog_ramp_speed = args.rot_ramp_speed)
    shop_bot_file.set_speed(z_move_speed = args.move_speed,
                            z_jog_speed = args.move_speed,
                            xy_move_speed = args.move_speed,
                            xy_jog_speed = args.move_speed,
                            b_move_speed = args.rot_move_speed,
                            b_jog_speed = args.rot_move_speed)
    
    
    points = [[0,0,0,0]]
    
    for i in range(20):
        for t in np.linspace(0, 1, 120):
            theta = t * t * 360.0
            x = 4*t
            y = 2*np.sqrt(t)
            z = 0.5*t
            scaled_theta = theta * 180.0 / 270.0
            points.append([x,y,z,scaled_theta])
        
    points.append([0,0,0,0])
    shop_bot_file.add_points(points)
    shop_bot_file.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    