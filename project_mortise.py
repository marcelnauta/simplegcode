# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 19:13:15 2020

@author: user
"""

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

    #parser.add_argument('--tenon-size-x', default=2.5, type=float,
    #                    help='Diameter of the holes being drilled.')
    
    #parser.add_argument('--tenon-size-y', default=1.5, type=float,
    #                    help='Diameter of the holes being drilled.')
    parser.add_argument('--tenon-size-x', default=4.25, type=float,
                        help='Diameter of the holes being drilled.')
    
    parser.add_argument('--tenon-size-y', default=1.0825, type=float,
                        help='Diameter of the holes being drilled.')
    
    parser.add_argument('--tenon-depth', default=2.0, type=float,
                        help='Diameter of the holes being drilled.')


    parser.add_argument('--ramp-speed', default=0.3, type=float,
                        help='Speed during corners in inches/second. Applies to tight radius helix only.')

    parser.add_argument('--move-speed', default=1.25, type=float,
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


    shop_bot_file = ShopBotFile('mortise.sbp')
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
    
    
    points = [[0,0,0]]
    x_width = args.tenon_size_x - args.bit_diameter
    y_width = args.tenon_size_y - args.bit_diameter
    scaling = 0.5
    x_passes = divide_into_equal_passes(-x_width/2.0, 0, scaling*args.bit_diameter, include_first = True)
    y_passes = divide_into_equal_passes(-y_width/2.0, 0, scaling*args.bit_diameter, include_first = True)
    z_passes = divide_with_clearout_depths(args.tenon_depth, 0.25, 3.0)
    print(x_passes)
    print(y_passes)
    print(z_passes)
    w = -args.bit_diameter*0.75/2
    for z in z_passes:
        for x, y in zip(x_passes[:1], y_passes[:1]):
            points.append([ x, y,-z])
            points.append([ x+w, y,-z])
            points.append([ x, y+w,-z])
            points.append([ x,-y-w,-z])
            points.append([ x+w,-y,-z])
            points.append([-x-w,-y,-z])
            points.append([-x,-y-w,-z])
            points.append([-x, y+w,-z])
            points.append([-x-w, y,-z])
            points.append([ x, y,-z])
        for x, y in zip(x_passes[1:], y_passes[1:]):
            points.append([ x, y,-z])
            points.append([ x,-y,-z])
            points.append([-x,-y,-z])
            points.append([-x, y,-z])
            points.append([ x, y,-z])
        
    points.append([0,0,0])
    shop_bot_file.add_points(points)
    shop_bot_file.close()
    
    points = np.array(points)
    width_x = 5
    width_y = 5
    origin_x = -0.5 * width_x
    origin_y = -0.5 * width_y
    workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.tenon_depth],
                                                   origin=[origin_x, origin_y, 0])

    time_per_hole = 20.0 # seconds (measured)
    path_estimate = workpiece_preview.estimate_time(points, args.move_speed)
    print('Estimating {:0.2f} minutes based on tool path'.format(path_estimate / 60.0))

    workpiece_preview.plot_three_axis(points, draw_bounding_box = False)  
    workpiece_preview.plot_birds_eye(points)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    