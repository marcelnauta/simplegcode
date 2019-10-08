
class RouterBit(object):
    def __init__(self, diameter_inches,
                       pass_depth_inches = 3/16.0, # Target cutting depth in a single pass
                       clearout_depth_inches = 1/2.0 + 1/32.0, # Maximum depth before a clearout operation
                       horizontal_pass_width_fraction = 0.4 # Target cutting width in a single pass as a fraction of bit diameter
                       ):
        self.diameter = diameter_inches
        self.radius = diameter_inches / 2.0
        self.pass_depth = pass_depth_inches
        self.clearout_depth = clearout_depth_inches
        self.pass_horiz = horizontal_pass_width_fraction * diameter_inches

    def get_radius(height_wrt_bottom):
        return self.radius