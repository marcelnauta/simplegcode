AXIS_X = 0
AXIS_Y = 1
AXIS_Z = 2
AXIS_ANGLE = 3

DIRECTION_CONVENTIONAL = 'conventional'
DIRECTION_CLIMB = 'climb'

SAFE_HEIGHT = 0.25

class MovementSpeeds(object):
    def __init__(self,
                 # Speeds are defined in inches/sec. (angles not specified in documentation)
                 xy_move_speed = 1.0,
                 z_move_speed = 1.0,
                 a_move_speed = 1.0,
                 b_move_speed = 1.0,
                 xy_jog_speed = 2.0,
                 z_jog_speed = 2.0,
                 a_jog_speed = 2.0,
                 b_jog_speed = 2.0,
                 ):
        self.xy_move_speed = xy_move_speed
        self.z_move_speed = z_move_speed
        self.a_move_speed = a_move_speed
        self.b_move_speed = b_move_speed

        self.xy_jog_speed = xy_jog_speed
        self.z_jog_speed = z_jog_speed
        self.a_jog_speed = a_jog_speed
        self.b_jog_speed = b_jog_speed

    def set(self,
            xy_move_speed = None,
            z_move_speed = None,
            a_move_speed = None,
            b_move_speed = None,
            xy_jog_speed = None,
            z_jog_speed = None,
            a_jog_speed = None,
            b_jog_speed = None,
            ):
        if xy_move_speed is not None:
            self.xy_move_speed = xy_move_speed
        if z_move_speed is not None:
            self.z_move_speed = z_move_speed
        if a_move_speed is not None:
            self.a_move_speed = a_move_speed
        if b_move_speed is not None:
            self.b_move_speed = b_move_speed

        if xy_jog_speed is not None:
            self.xy_jog_speed = xy_jog_speed
        if z_jog_speed is not None:
            self.z_jog_speed = z_jog_speed
        if a_jog_speed is not None:
            self.a_jog_speed = a_jog_speed
        if b_jog_speed is not None:
            self.b_jog_speed = b_jog_speed

    def get_command(self):
        values = [self.xy_move_speed,
                  self.z_move_speed,
                  self.a_move_speed,
                  self.b_move_speed,

                  self.xy_jog_speed,
                  self.z_jog_speed,
                  self.a_jog_speed,
                  self.b_jog_speed,
                  ]
        return ','. join(['VS'] + ['{0:.2f}'.format(v) for v in values]) + '\n'

class RampSpeeds(object):
    def __init__(self,
                 # The Move Ramp Speed sets the start and end speed for acceleration
                 # and deceleration ramps for each axis for cutting speeds.
                 # In other words, the slowest move speed in each axis.
                 xy_move_ramp_speed = 0.4,
                 z_move_ramp_speed = 0.4,
                 a_move_ramp_speed = 0.4,
                 b_move_ramp_speed = 0.4,
                 # The Jog Ramp Speed sets the start and end speed for acceleration
                 # and deceleration ramps for each axis for rapid transit speeds.
                 # In other words, the slowest jog speed in each axis.
                 xy_jog_ramp_speed = 0.4,
                 z_jog_ramp_speed = 0.4,
                 a_jog_ramp_speed = 0.4,
                 b_jog_ramp_speed = 0.4,
                 # The distances over which a change of 2 Units is made in Move Speed
                 # or Jog speed during ramping. This parameter thus sets the
                 # acceleration rate and deceleration rate.
                 move_ramp_rate = 0.2,
                 # This value adjusts the sensitivity of 3D ramping. Decreasing makes
                 # 3D ramping more sensitive, increasing makes it less sensitive.
                 threshold_3d = 100.0,
                 # The minimum distance considered in making ramping decisions.
                 minimum_distance = 0.15,
                 #The percentage of speed reduction when ramping is triggered for a
                 # slight change in direction.Angular changes in direction produce
                 # increasingly greater reduction in speed from this value to 100%
                 # of the full ramp reduction.
                 slow_corner_speed = 65.0,
                 # Sets the acceleration rate for KeyPad mode
                 keypad_ramp_rate = 0.2,
                 # Circle diameter below which execution will occur at Ramp Speed
                 # rather thanMove Speed. Small circles would be executed too
                 # abruptly if they accelerated to full Move Speed within their
                 # small diameter. This would be a particular issue with a spiral
                 # plunge.
                 small_circle_diameter = 0.25,
                 ):
        self.xy_move_ramp_speed = xy_move_ramp_speed
        self.z_move_ramp_speed = z_move_ramp_speed
        self.a_move_ramp_speed = a_move_ramp_speed
        self.b_move_ramp_speed = b_move_ramp_speed

        self.xy_jog_ramp_speed = xy_jog_ramp_speed
        self.z_jog_ramp_speed = z_jog_ramp_speed
        self.a_jog_ramp_speed = a_jog_ramp_speed
        self.b_jog_ramp_speed = b_jog_ramp_speed

        self.move_ramp_rate = move_ramp_rate
        self.threshold_3d = threshold_3d
        self.minimum_distance = minimum_distance
        self.slow_corner_speed = slow_corner_speed
        self.keypad_ramp_rate = keypad_ramp_rate
        self.small_circle_diameter = small_circle_diameter

    def set(self,
            xy_move_ramp_speed = None,
            z_move_ramp_speed = None,
            a_move_ramp_speed = None,
            b_move_ramp_speed = None,
            xy_jog_ramp_speed = None,
            z_jog_ramp_speed = None,
            a_jog_ramp_speed = None,
            b_jog_ramp_speed = None,
            move_ramp_rate = None,
            threshold_3d = None,
            minimum_distance = None,
            slow_corner_speed = None,
            keypad_ramp_rate = None,
            small_circle_diameter = None,
            ):
        if xy_move_ramp_speed is not None:
            self.xy_move_ramp_speed = xy_move_ramp_speed
        if z_move_ramp_speed is not None:
            self.z_move_ramp_speed = z_move_ramp_speed
        if a_move_ramp_speed is not None:
            self.a_move_ramp_speed = a_move_ramp_speed
        if b_move_ramp_speed is not None:
            self.b_move_ramp_speed = b_move_ramp_speed

        if xy_jog_ramp_speed is not None:
            self.xy_jog_ramp_speed = xy_jog_ramp_speed
        if z_jog_ramp_speed is not None:
            self.z_jog_ramp_speed = z_jog_ramp_speed
        if a_jog_ramp_speed is not None:
            self.a_jog_ramp_speed = a_jog_ramp_speed
        if b_jog_ramp_speed is not None:
            self.b_jog_ramp_speed = b_jog_ramp_speed

        if move_ramp_rate is not None:
            self.move_ramp_rate = move_ramp_rate
        if threshold_3d is not None:
            self.threshold_3d = threshold_3d
        if minimum_distance is not None:
            self.minimum_distance = minimum_distance
        if slow_corner_speed is not None:
            self.slow_corner_speed = slow_corner_speed
        if keypad_ramp_rate is not None:
            self.keypad_ramp_rate = keypad_ramp_rate
        if small_circle_diameter is not None:
            self.small_circle_diameter = small_circle_diameter

    def get_command(self):
        values = [self.xy_move_ramp_speed,
                  self.z_move_ramp_speed,
                  self.a_move_ramp_speed,
                  self.b_move_ramp_speed,

                  self.xy_jog_ramp_speed,
                  self.z_jog_ramp_speed,
                  self.a_jog_ramp_speed,
                  self.b_jog_ramp_speed,

                  self.move_ramp_rate,
                  self.threshold_3d,
                  self.minimum_distance,
                  self.slow_corner_speed,
                  # Manual is inconsistent with these parameters
                  #self.keypad_ramp_rate,
                  #self.small_circle_diameter,
                  ]
        return ','. join(['VR'] + ['{0:.2f}'.format(v) for v in values]) + '\n'
