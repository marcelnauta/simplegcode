# Finger joint generator for Rhinocerus CAD
# Copyright (C) 2021 Marcel Nauta
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from collections import OrderedDict
import copy
import json
import math
import os
import sys

# You likely need your own install path, but this might help find it
try:
    rhino_cwd = os.getcwd()
    rhino_plugins_path = os.path.join(rhino_cwd, '..', 'Plug-ins')
    for subdir in os.listdir(rhino_plugins_path):
        if subdir.startswith('IronPython'):
            sys.path.append(os.path.join(rhino_plugins_path, subdir, 'settings', 'lib', 'rhinoscript'))
except:
    pass

try:
    import rhinoscriptsyntax as rs
    import Rhino
    import utility as rhutil
    import scriptcontext
except ImportError:
    rs = None
    Rhino = None
    rhutil = None
    scriptcontext = None


DEFAULT_BOARD_THICKNESS = 4 # mm
DEFAULT_BOARD_WIDTH = 150 # mm
DEFAULT_BOARD_LENGTH = 300 # mm

# Helper functions for building deltas. Mostly matches functionality from numpy,
# but numpy is not shipped with the Python in Rhino
def linspace(min=0.0, max=1.0, num=51):
    delta = (max-min) / (num-1.0)
    return [min + i*delta for i in range(num)]

class Extremum(object):
    def __init__(self, axis, minmax):
        self.axis = axis
        self.minmax = minmax

    def is_min(self):
        return self.minmax == 'min'

    def to_single(self):
        return '{}_{}'.format(self.axis, self.minmax)

    def to_pair(self):
        return (self.axis, self.minmax)

    def to_axis_index(self):
        return ALL_AXES.index(self.axis)

    def to_minmax_sign(self):
        return -1 if self.is_min() else 1

    def perpendicular_axis(self, other_extremum):
        mapping = {
            'xy': 'z',
            'yx': 'z',
            'xz': 'y',
            'zx': 'y',
            'zy': 'x',
            'yz': 'x',
            }
        return mapping[self.axis + other_extremum.axis]

    def opposite(self):
        return Extremum(self.axis, 'max' if self.is_min() else 'min')

    def __hash__(self):
        return hash(self.to_single())

    def __str__(self):
        return self.to_single()

    def __repr__(self):
        return self.to_single()

    def __eq__(self, other):
        if isinstance(other, Extremum):
            return self.to_single() == other.to_single()
        elif isinstance(other, str):
            return self.to_single() == other
        return False

    def __lt__(self, other):
        if isinstance(other, Extremum):
            return self.to_single() < other.to_single()
        elif isinstance(other, str):
            return self.to_single() < other
        return False

ALL_AXES = ['x', 'y', 'z']
ALL_MINMAX = ['min', 'max']
ALL_EXTREMA = []
for axis in ALL_AXES:
    for minmax in ALL_MINMAX:
        ALL_EXTREMA.append(Extremum(axis, minmax))

def dict_to_point(point_as_dict):
    return [point_as_dict[axis] for axis in ALL_AXES]

def get_other_axes(axis):
    return [other_axis for other_axis in ALL_AXES if not other_axis == axis]

# Rhino draws boxes with an annoying syntax
class RhinoBox(object):
    def get_midpoint(self):
        midpoint = []
        for axis in ALL_AXES:
            axis_midpoint = 0
            for minmax in ALL_MINMAX:
                axis_midpoint += 0.5 * self.extreme_points[Extremum(axis, minmax)]
            midpoint.append(axis_midpoint)
        return midpoint


    def __init__(self, extreme_points):
        self.extreme_points = extreme_points
        self.x_min = extreme_points['x_min']
        self.x_max = extreme_points['x_max']
        self.y_min = extreme_points['y_min']
        self.y_max = extreme_points['y_max']
        self.z_min = extreme_points['z_min']
        self.z_max = extreme_points['z_max']

    def draw(self):
        points = []
        for z_point in [self.z_min, self.z_max]:
          for x_point, y_point in [[self.x_max, self.y_max], [self.x_min, self.y_max], [self.x_min, self.y_min], [self.x_max, self.y_min]]:
              points.append((x_point, y_point, z_point))
        if rhutil is not None:
            corners = points
            # Copied from surface.py and remove the redraw to save time
            box = rhutil.coerce3dpointlist(corners, True)
            brep = Rhino.Geometry.Brep.CreateFromBox(box)
            if not brep: raise ValueError("unable to create brep from box")
            rc = scriptcontext.doc.Objects.AddBrep(brep)
            return rc
        else:
            return rs.AddBox(points)

class Rotation(object):
    def __init__(self, enabled=False, normal_direction=[0,0,0], angle_degrees=90):
        self.enabled = enabled
        self.normal_direction = normal_direction
        if normal_direction == 'x':
            self.normal_direction = [1.0, 0.0, 0.0]
        elif normal_direction == 'y':
            self.normal_direction = [0.0, 1.0, 0.0]
        elif normal_direction == 'z':
            self.normal_direction = [0.0, 0.0, 1.0]
        self.angle_degrees = angle_degrees

    def rotate_objects(self, object_ids, center_point):
       if not self.enabled:
           return object_ids
       return rs.RotateObjects(object_ids, center_point, self.angle_degrees, self.normal_direction, copy=False)

    def rotate_objects_back(self, object_ids, center_point):
       if not self.enabled:
           return object_ids
       return rs.RotateObjects(object_ids, center_point, -self.angle_degrees, self.normal_direction, copy=False)

class BasicBoard(object):
    def set_rotations_by_axes(self, normal, grain, is_flipped):
        if normal == grain:
            raise ValueError("Normal axis and grain axis cannot point in the same direction")

        if is_flipped:
            self.rotations.append(Rotation(enabled=True, normal_direction='y', angle_degrees = 180))

        if normal == 'z':
            # No need to rotate the normal, only rotate if grain is along y
            if grain == 'y':
                self.rotations.append(Rotation(enabled=True, normal_direction='z'))
        elif normal == 'y':
            # To rotate z to y, we need a rotation about x
            self.rotations.append(Rotation(enabled=True, normal_direction='x'))
            # y and z were swapped, therefore grain is currently along z
            if grain == 'y':
                self.rotations.append(Rotation(enabled=True, normal_direction='y'))
        elif normal == 'x':
            # To rotate z to x, we need a rotation about y
            self.rotations.append(Rotation(enabled=True, normal_direction='y'))
            # x and z were swapped, therefore grain is currently along z
            if grain == 'y':
                self.rotations.append(Rotation(enabled=True, normal_direction='x'))
        else:
            ValueError('Unrecognized normal axis. Recieved {}'.format(normal))

    def __init__(self, thickness, width, length, normal, grain, is_flipped, center, num_boards_wide = 1, num_boards_long = 2):
        self.rhino_boxes = []
        xlines = linspace(-0.5 * length, 0.5 * length, num_boards_long+1)
        ylines = linspace(-0.5 * width, 0.5 * width, num_boards_wide+1)
        overlap = thickness
        for iLength in range(num_boards_long):
            for iWidth in range(num_boards_wide):
                extreme_points = OrderedDict()
                extreme_points['x_min'] = xlines[iLength]   - (0.5 * overlap if iLength != 0                 else 0.0)
                extreme_points['x_max'] = xlines[iLength+1] + (0.5 * overlap if iLength != num_boards_long-1 else 0.0)
                extreme_points['y_min'] = ylines[iWidth]    - (0.5 * overlap if iWidth != 0                  else 0.0)
                extreme_points['y_max'] = ylines[iWidth+1]  + (0.5 * overlap if iWidth != num_boards_wide-1  else 0.0)
                extreme_points['z_min'] = -0.5 * thickness
                extreme_points['z_max'] =  0.5 * thickness
                self.rhino_boxes.append(RhinoBox(extreme_points))
        self.rotations = []
        self.set_rotations_by_axes(normal, grain, is_flipped)
        self.thickness = thickness
        self.width = width
        self.length = length
        self.center = center
        self.num_boards_wide = num_boards_wide
        self.num_boards_long = num_boards_long
        self.overlap = overlap
        self.xlines = xlines
        self.ylines = ylines

    def draw(self, width_finger_generator, length_finger_generator):
        object_ids = [rhino_box.draw() for rhino_box in self.rhino_boxes]

        locations = length_finger_generator.get_locations()
        for iLength in range(self.num_boards_long):
            for iWidth in range(self.num_boards_wide-1):
                joint_ids = []
                if length_finger_generator.min_inset > 0:
                    joint_extremes = OrderedDict()
                    joint_extremes['y_min'] = self.ylines[iWidth+1]
                    joint_extremes['y_max'] = self.ylines[iWidth+1] + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('x', 'min')] = -0.5*self.length
                    joint_extremes[Extremum('x', 'max')] = locations[0] - 0.5*self.length
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                if length_finger_generator.max_inset > 0:
                    joint_extremes = OrderedDict()
                    joint_extremes['y_min'] = self.ylines[iWidth+1]
                    joint_extremes['y_max'] = self.ylines[iWidth+1] + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('x', 'min')] = locations[-1] - 0.5*self.length
                    joint_extremes[Extremum('x', 'max')] = 0.5*self.length
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                for i in range(1,len(locations),2):
                    joint_extremes = OrderedDict()
                    joint_extremes['y_min'] = self.ylines[iWidth+1]  - 0.5 * self.overlap
                    joint_extremes['y_max'] = self.ylines[iWidth+1]  + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('x', 'min')] = locations[i-1] - 0.5*self.length
                    joint_extremes[Extremum('x', 'max')] = locations[i  ] - 0.5*self.length
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                primary_board_id = iLength * self.num_boards_wide + iWidth
                secondary_board_id = iLength * self.num_boards_wide + (iWidth+1)
                object_ids[primary_board_id] = rs.BooleanDifference([object_ids[primary_board_id]], joint_ids, True)[0]
                temp_board_id = rs.BooleanDifference([object_ids[secondary_board_id]], [object_ids[primary_board_id]], False)[0]
                rs.DeleteObject(object_ids[secondary_board_id])
                object_ids[secondary_board_id] = temp_board_id

        locations = width_finger_generator.get_locations()
        for iLength in range(self.num_boards_long-1):
            for iWidth in range(self.num_boards_wide):
                joint_ids = []
                if width_finger_generator.min_inset > 0:
                    joint_extremes = OrderedDict()
                    joint_extremes['x_min'] = self.xlines[iLength+1]
                    joint_extremes['x_max'] = self.xlines[iLength+1] + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('y', 'min')] = -0.5*self.width
                    joint_extremes[Extremum('y', 'max')] = locations[0] - 0.5*self.width
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                if width_finger_generator.max_inset > 0:
                    joint_extremes = OrderedDict()
                    joint_extremes['x_min'] = self.xlines[iLength+1]
                    joint_extremes['x_max'] = self.xlines[iLength+1] + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('y', 'min')] = locations[-1] - 0.5*self.width
                    joint_extremes[Extremum('y', 'max')] = 0.5*self.width
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                for i in range(1,len(locations),2):
                    joint_extremes = OrderedDict()
                    joint_extremes['x_min'] = self.xlines[iLength+1]  - 0.5 * self.overlap
                    joint_extremes['x_max'] = self.xlines[iLength+1]  + 0.5 * self.overlap
                    joint_extremes['z_min'] = -0.5 * self.thickness
                    joint_extremes['z_max'] =  0.5 * self.thickness
                    joint_extremes[Extremum('y', 'min')] = locations[i-1] - 0.5*self.width
                    joint_extremes[Extremum('y', 'max')] = locations[i  ] - 0.5*self.width
                    joint_box = RhinoBox(joint_extremes)
                    joint_ids.append(joint_box.draw())
                primary_board_id = iLength * self.num_boards_wide + iWidth
                secondary_board_id = (iLength+1) * self.num_boards_wide + iWidth
                object_ids[primary_board_id] = rs.BooleanDifference([object_ids[primary_board_id]], joint_ids, True)[0]
                temp_board_id = rs.BooleanDifference([object_ids[secondary_board_id]], [object_ids[primary_board_id]], False)[0]
                rs.DeleteObject(object_ids[secondary_board_id])
                object_ids[secondary_board_id] = temp_board_id

        rotation_center = [0.0, 0.0, 0.0]
        for rotation in self.rotations:
            object_ids = rotation.rotate_objects(object_ids, rotation_center)
        translation = []
        for i, axis in enumerate(ALL_AXES):
            translation.append(self.center[axis] - rotation_center[i])
        if translation != [0.0, 0.0, 0.0]:
            object_ids = rs.MoveObjects(object_ids, translation)
        return object_ids

    def laydown(self, object_ids):
        rotation_center = [0.0, 0.0, 0.0]
        translation = []
        for i, axis in enumerate(ALL_AXES):
            translation.append(- self.center[axis] + rotation_center[i])
        if translation != [0.0, 0.0, 0.0]:
            object_ids = rs.MoveObjects(object_ids, translation)
        for rotation in self.rotations[::-1]:
            object_ids = rotation.rotate_objects_back(object_ids, rotation_center)
        if 0:
            # No known way to get number of surfaces and know which is which...
            for i, object_id in enumerate(object_ids):
                brep = rhutil.coercebrep(object_id, True)
                for index in range(brep.Faces.Count):
                    face = brep.Faces[index]
                    newbrep = face.DuplicateFace(True)
                    surface_id = scriptcontext.doc.Objects.AddBrep(newbrep)
                    surface = rhutil.coercesurface(surface_id, True)
                    surface_normal = surface.NormalAt(0.5, 0.5)
                    if surface_normal[2] > 0.99: # positive surface normal in z direction (should only find one)
                        break
                    else:
                        scriptcontext.doc.Objects.Delete(surface_id, True)
                rs.DeleteObject(object_id)
                object_ids[i] = rhutil.coerceguid(surface_id)
            # END of magic, each object_id should now hold a single surface
            for i, object_id in enumerate(object_ids):
                line_ids = rs.DuplicateEdgeCurves(object_id)
                rs.DeleteObject(object_id)
                object_ids[i] = rs.JoinCurves(line_ids, delete_input=True)
        return object_ids

################################################################################

class ConfigEntry(object):
    def __init__(self, key, description = '', default = None):
        self.key = key
        self.description = description
        self.default = default
        self.value = None

    def from_property_tree(self, ptree):
        try:
            self.value = ptree[self.key]
        except KeyError:
            print('Unable to find {} for {} in {}'.format(self.key, self.description, ptree))
            raise
        return self.value

    def from_user_prompt(self):
        return None

    def to_property_tree(self, ptree):
        if self.value is None:
            raise NotImplementedError('self.value must be set or a to_property_tree() overridden')
        ptree[self.key] = self.value

    def get(self):
        if self.value is None:
            raise NotImplementedError('self.value must be set or a get() overridden')
        return self.value

    def __str__(self):
        individual_ptree = {}
        self.to_property_tree(individual_ptree)
        return str(individual_ptree)


class FloatEntry(ConfigEntry):
    def __init__(self, key, description, default, minimum = None, maximum = None):
        ConfigEntry.__init__(self, key, description, default)
        self.minimum = minimum
        self.maximum = maximum

    def from_user_prompt(self):
        self.value = rs.GetReal(message=self.description, number=self.default, minimum=self.minimum, maximum=self.maximum)
        return self.value

class IntegerEntry(ConfigEntry):
    def __init__(self, key, description, default, minimum = None, maximum = None):
        ConfigEntry.__init__(self, key, description, default)
        self.minimum = minimum
        self.maximum = maximum

    def from_user_prompt(self):
        self.value = rs.GetInteger(message=self.description, number=self.default, minimum=self.minimum, maximum=self.maximum)
        return self.value

class BoolEntry(ConfigEntry):
    def from_user_prompt(self):
        rhino_rv = rs.GetBoolean(message=self.description, items=[('Enabled','no','yes')], defaults=[self.default])
        self.value = rhino_rv[0] if rhino_rv is not None else None
        return self.value

class AxisEntry(ConfigEntry):
    def from_user_prompt(self):
        self.value = rs.GetString(message=self.description, defaultString=self.default,
                                  strings=ALL_AXES)
        return self.value

class EntryFactory(object):
    def __init__(self, factory_entry, factory_fn, factory_kwargs = {}):
        self.factory_entry = factory_entry
        self.factory_fn = factory_fn
        self.factory_kwargs = factory_kwargs

    def create_entry(self, class_type, key):
        self.derived_entry = self.factory_fn(class_type, key, **self.factory_kwargs)
        return self.derived_entry

    def get_entry(self):
        return self.derived_entry

class ParentEntry(ConfigEntry):
    def __init__(self, *args, **kwargs):
        ConfigEntry.__init__(self, *args, **kwargs)
        self.children = OrderedDict()
        self.factories = OrderedDict()

    def _children_to_attributes(self):
        for child_key in self.children:
            setattr(self, str(child_key), self.children[child_key].value)
        #for factory_key in self.factories:
        #    setattr(self, str(factory_key), self.factories[factory_key].value)

    def from_user_prompt(self):
        self.value = OrderedDict()
        for key in self.children:
            self.value[key] = self.children[key].from_user_prompt()
        for key in self.factories:
            factory_entry = self.factories[key].factory_entry
            # The factory entry must return a single value such as a string or integer
            class_type = factory_entry.from_user_prompt() if factory_entry is not None else None
            derived_class_entry = self.factories[key].create_entry(class_type, key)
            # The derived class must return a dictionary
            factory_value = derived_class_entry.from_user_prompt()
            factory_value[factory_entry.key] = class_type
            self.value[key] = factory_value
        self._children_to_attributes()
        return self.value

    def from_property_tree(self, ptree):
        self.value = OrderedDict()
        for chlid_key in self.children:
            self.value[chlid_key] = self.children[chlid_key].from_property_tree(ptree[self.key])
        for key in self.factories:
            factory_entry = self.factories[key].factory_entry
            # The factory entry must return a single value such as a string or integer
            class_type = factory_entry.from_property_tree(ptree[self.key][key]) if factory_entry is not None else None
            derived_class_entry = self.factories[key].create_entry(class_type, key)
            # The derived class must return a dictionary
            factory_value = derived_class_entry.from_property_tree(ptree[self.key])
            factory_value[factory_entry.key] = class_type
            self.value[key] = factory_value
        self._children_to_attributes()
        return self.value

class PointEntry(ParentEntry):
    def __init__(self, *args, **kwargs):
        ParentEntry.__init__(self, *args, **kwargs)
        for axis in ALL_AXES:
            self.children[axis] = FloatEntry(axis, description = '{} along {} axis'.format(self.description, axis),
                                             default = self.default if self.default is not None else 0.0)
    def to_point(self):
        return [self.children[axis].value for axis in ALL_AXES]

def point_factory(class_type, key, **kwargs):
    return PointEntry(key)

class PointFactory(EntryFactory):
    def __init__(self, factory_entry, **kwargs):
        EntryFactory.__init__(self, factory_entry = factory_entry,
                              factory_fn = point_factory,
                              factory_kwargs = kwargs)

##################### FINGERS ######################

class FingerTypeEntry(ConfigEntry):
    def from_user_prompt(self):
        self.value = rs.GetString(message='Finger type for ' + self.description, defaultString='New',
                                  strings=['Custom', 'SquishIn', 'SquishOut', 'SquishMin', 'SquishMax', 'SquishMulti', 'Equal', 'Random', 'New'])
        return self.value

class FingerLocations(ParentEntry):
    def __init__(self, *args, **kwargs):
        points_from_json = False
        if 'points_from_json' in kwargs:
            points_from_json = kwargs.pop('points_from_json')
        ParentEntry.__init__(self, *args, **kwargs)
        self.children['min_size'] = FloatEntry('min_size', description = 'Min finger size', default = 2.0)
        self.children['min_inset'] = FloatEntry('min_inset', description = 'Min inset (likely thickness)', default = 4.0)
        self.children['max_inset'] = FloatEntry('max_inset', description = 'Max inset (likely thickness)', default = 4.0)
        if points_from_json:
            self.children['min_point'] = PointEntry('min_point', description = 'Start point', default = 0.0)
            self.children['max_point'] = PointEntry('max_point', description = 'End point', default = 1.0)
        self.curve_id = None
        self.intersector_size = 1000.0

    def set_curve(self, curve_id):
        self.curve_id = curve_id
        # Set the "length" to be the return value of Rhino's curve domain. Appears to be in distance units
        curve_domain = rs.CurveDomain(curve_id)
        self.length = curve_domain[1] - curve_domain[0]

    def get_total_length(self):
        return self.length - self.min_inset - self.max_inset

    def init_with_user_curve(self):
        self.set_curve(rs.GetCurveObject('Select a curve to base fingers on')[0])

    def get_locations(self, draw_type = None):
        locations = self._get_locations()
        if len(locations) == 0:
            raise ValueError('No locations set')
        if sorted(locations) != locations:
            raise ValueError('Locations are not in ascending order')
        if locations[0] <= self.min_inset:
            raise ValueError('First point is outside the bounds of min inset')
        if locations[-1] >= (self.length-self.max_inset):
            raise ValueError('First point is outside the bounds of max inset')
        if draw_type is not None:
            new_locations = []
            if draw_type == 'mini':
                for i in range(1,len(locations)-1):
                    direction = +1 if i % 2 == 0 else -1
                    i_shift = i + direction
                    if i_shift > 0 and i_shift < len(locations):
                        new_locations.append(2*locations[i]/3 + 1*locations[i_shift]/3)
            if draw_type == 'mini_offset':
                for i in range(len(locations)):
                    direction = +1 if i % 2 == 1 else -1
                    i_shift = i + direction
                    if i_shift > 0 and i_shift < len(locations):
                        new_locations.append(2*locations[i]/3 + 1*locations[i_shift]/3)
            locations = new_locations
        return locations


    def _draw_plane(self, curve_id, curve_t, plane_size):
        point_on_curve = rs.EvaluateCurve(curve_id, curve_t)
        surface_normal = rs.CurveTangent(curve_id, curve_t)
        plane_info = rs.PlaneFromNormal(point_on_curve, surface_normal)
        surface_id = rs.AddPlaneSurface(plane_info, plane_size, plane_size)
        # The plane is drawn starting from the bottom right, so move it to start at the center
        surface_domain_midpoint = []
        for domain_dir in [0, 1]:
            surface_domain = rs.SurfaceDomain(surface_id, domain_dir)
            surface_domain_midpoint.append(0.5 * (surface_domain[0] + surface_domain[1]))
        surface_midpoint = rs.EvaluateSurface(surface_id, surface_domain_midpoint[0], surface_domain_midpoint[1])
        surface_translation = point_on_curve - surface_midpoint
        surface_id = rs.MoveObject(surface_id, surface_translation)
        return surface_id

    def _make_curve_segment(self, curve_id, curve_t_start, curve_t_end):
        first_point_on_curve = rs.EvaluateCurve(curve_id, curve_t_start)
        last_point_on_curve = rs.EvaluateCurve(curve_id, curve_t_end)
        return rs.AddCurve([first_point_on_curve, last_point_on_curve])

    def _extrude_along_segment(self, curve_id, curve_t_start, curve_t_end, intersect_width):
        curve_segment_id = self._make_curve_segment(curve_id, curve_t_start, curve_t_end)
        plane_id = self._draw_plane(curve_segment_id, rs.CurveDomain(curve_segment_id)[0], intersect_width)
        new_brep_id = rs.ExtrudeSurface(plane_id, curve_segment_id)
        rs.DeleteObjects([curve_segment_id, plane_id])
        return new_brep_id

    def get_min_point(self):
        return self.children['min_point'].to_point()

    def get_max_point(self):
        return self.children['max_point'].to_point()

    def draw(self, draw_type = None):
        INTERSECT_WIDTH = self.intersector_size
        # Create a polyline
        curve_id = self.curve_id
        curve_domain = rs.CurveDomain(curve_id)
        locations = self.get_locations(draw_type)
        finger_brep_ids = []
        for iLocation in range(1, len(locations), 2):
            curve_t_neg = locations[iLocation-1] + curve_domain[0]
            curve_t_pos = locations[iLocation  ] + curve_domain[0]
            finger_brep_ids.append(self._extrude_along_segment(curve_id, curve_t_neg, curve_t_pos, INTERSECT_WIDTH))
        return finger_brep_ids

    def draw_old(self, draw_type = None):
        INTERSECT_WIDTH = self.intersector_size
        # Create a polyline
        curve_id = self.curve_id
        curve_domain = rs.CurveDomain(curve_id)
        locations = self.get_locations(draw_type)
        start_plane_id = self._draw_plane(curve_id, curve_domain[0], INTERSECT_WIDTH)
        last_test_plane_id = self._draw_plane(curve_id, 0.5*(locations[-1]+curve_domain[0] + curve_domain[1]), INTERSECT_WIDTH)
        main_brep_id = rs.ExtrudeSurface(start_plane_id, curve_id)
        #main_brep_id = rs.AddPipe(curve_id, curve_domain, [0.5*INTERSECT_WIDTH, 0.5*INTERSECT_WIDTH], cap=1)
        if main_brep_id is None:
            raise RuntimeError('Failed to make initial brep for fingers')
        finger_brep_ids = []
        removable_ids = [start_plane_id, last_test_plane_id, main_brep_id]
        for iLocation in range(1, len(locations), 2):
            curve_t_neg = locations[iLocation-1] + curve_domain[0]
            curve_t_pos = locations[iLocation  ] + curve_domain[0]
            curve_t_mid = 0.5*(locations[iLocation] + locations[iLocation-1]) + curve_domain[0]
            point_on_curve = rs.EvaluateCurve(curve_id, curve_t_mid)
            cut_plane_id_neg = self._draw_plane(curve_id, curve_t_neg, 1*INTERSECT_WIDTH)
            cut_plane_id_pos = self._draw_plane(curve_id, curve_t_pos, 1*INTERSECT_WIDTH)
            neg_brep_ids = rs.SplitBrep(main_brep_id, cut_plane_id_neg, delete_input=False)
            removable_ids.extend(neg_brep_ids)
            desirable_id = None
            # Not sure how to find the brep we want except to cut with another plane and test
            # if the intersection was None
            for neg_brep_id in neg_brep_ids:
                pos_brep_ids = rs.SplitBrep(neg_brep_id, cut_plane_id_pos, delete_input=False)
                if pos_brep_ids is not None:
                    for pos_brep_id in pos_brep_ids:
                        test_brep_ids = rs.SplitBrep(pos_brep_id, last_test_plane_id, delete_input=False)
                        if test_brep_ids is not None:
                            removable_ids.append(pos_brep_id)
                            removable_ids.extend(test_brep_ids)
                        else:
                            desirable_id = pos_brep_id
            if desirable_id is None:
                raise RuntimeError('Unable to create a brep for finger')
            finger_brep_id = rs.JoinSurfaces([desirable_id, cut_plane_id_neg, cut_plane_id_pos],  delete_input=True)
            if finger_brep_id is not None:
                finger_brep_ids.append(finger_brep_id)
            else:
                raise RuntimeError('Unable to join surfaces')
        rs.DeleteObjects(removable_ids)
        return finger_brep_ids

    def apply(self, primary_object_ids, secondary_object_ids):
        finger_brep_ids = self.draw()
        overlap_ids = rs.BooleanIntersection(primary_object_ids, secondary_object_ids, delete_input=False)
        overlap_ids = rs.BooleanIntersection(overlap_ids, finger_brep_ids, delete_input=True)
        secondary_object_ids = rs.BooleanDifference(secondary_object_ids, overlap_ids, delete_input=True)
        old_primary_object_ids = primary_object_ids
        primary_object_ids = rs.BooleanDifference(primary_object_ids, secondary_object_ids, delete_input=False)
        rs.DeleteObjects(old_primary_object_ids)
        return primary_object_ids, secondary_object_ids

    def apply_in_groups(self, primary_object_ids_dict, secondary_object_ids_dict, draw_type = None, polarity = 0):
        '''Changes the input dictionaries and fills in the new object ids.'''
        if len(primary_object_ids_dict) == 0:
            return None
        if len(secondary_object_ids_dict) == 0:
            return None
        if polarity:
            tmp_dict = primary_object_ids_dict
            primary_object_ids_dict = secondary_object_ids_dict
            secondary_object_ids_dict = tmp_dict
        # Make all of the primary objects from the fingers
        finger_brep_ids = self.draw(draw_type)
        select_object_ids = [x for sublist in primary_object_ids_dict.values() for x in sublist]
        all_secondary_object_ids = [x for sublist in secondary_object_ids_dict.values() for x in sublist]
        select_object_ids.extend(all_secondary_object_ids)
        #rs.GetObjects('pause', select=True, objects=select_object_ids)
        for primary_key in primary_object_ids_dict:
            primary_object_ids = primary_object_ids_dict[primary_key]
            overlap_ids = rs.BooleanIntersection(primary_object_ids, all_secondary_object_ids, delete_input=False)
            finger_ids = rs.BooleanIntersection(overlap_ids, finger_brep_ids, delete_input=False)
            rs.DeleteObjects(overlap_ids)
            primary_object_ids = rs.BooleanDifference(primary_object_ids, finger_ids, delete_input=True)
            primary_object_ids_dict[primary_key] = primary_object_ids
        rs.DeleteObjects(finger_brep_ids)
        # Make all of the secondary objects from the promary ones
        all_primary_object_ids = [x for sublist in primary_object_ids_dict.values() for x in sublist]
        for secondary_key in secondary_object_ids_dict:
            secondary_object_ids = secondary_object_ids_dict[secondary_key]
            # My understanding is that we should be able to use one call to BooleanDifference, but
            # in some cases, if one of the objects does not intersect, then the function returns None,
            # even though other objects have a valid intersection. Therefore, we loop (tested in Rhino 7)
            new_ids = []
            for secondary_object_id in secondary_object_ids:
                tmp_ids = rs.BooleanDifference([secondary_object_id], all_primary_object_ids, delete_input=False)
                if tmp_ids is not None:
                    new_ids.extend(tmp_ids)
            if len(new_ids) == 0:
                print secondary_object_ids, all_primary_object_ids
                rs.ObjectColor(secondary_object_ids, [0,255,0])
                rs.ObjectColor(all_primary_object_ids, [0,0,255])
                raise RuntimeError('Unable to make secondary objects for figners')
            else:
                rs.DeleteObjects(secondary_object_ids)
                secondary_object_ids_dict[secondary_key] = new_ids

    def apply_user(self):
        self.intersector_size = rs.GetReal('Radius of intersection pipe', number=1000.0, minimum=0)
        primary_object_ids = rs.GetObjects('Select primary objects', filter=16)
        secondary_object_ids = rs.GetObjects('Select secondary objects', filter=16)
        self.apply(primary_object_ids, secondary_object_ids)


#class CustomFingerLocations(FingerLocations):
#    def __init__(self, *args, **kwargs):
#        FingerLocations.__init__(self, *args, **kwargs)
#        min_default = 0.25 * (self.length - self.max_inset - self.min_inset) + self.min_inset
#        max_default = 0.75 * (self.length - self.max_inset - self.min_inset) + self.min_inset
#        far = self.min_size * 1.5
#        near = self.min_size * 0.5
#        defaults = [min_default-far, min_default-near, min_default+near, min_default+far,
#                    max_default-far, max_default-near, max_default+near, max_default+far]
#
#        self.locations_str = rs.GetString(message='Joint transition points in mm',
#                                          defaultString=','.join(str(x) for x in defaults))
#
#    def _get_locations(self):
#        return [float(x.strip()) for x in self.locations_str.split(',')]

class VariableFingerLocations(FingerLocations):
    def __init__(self, *args, **kwargs):
        FingerLocations.__init__(self, *args, **kwargs)
        self.children['max_size'] = FloatEntry('max_size', description = 'Max finger size', default = 6.0)

    def get_scalar(self, normalized_location):
        ''' Return a value between 0.0 and 1.0 that represents finger width for the normalized_location. '''
        return 0.0

    def get_delta(self, normalized_location):
        return self.min_size + self.get_scalar(normalized_location) * (self.max_size-self.min_size)

    def re_seed(self):
        ''' Things like random number generators need to be re-seeded every iteration. '''
        pass

    def _get_locations(self):
        # Just because I don't trust Rhino's checks and might take these numbers in another way someday
        if self.min_size <= 0.0:
            raise ValueError('min_size must be strictly positive')
        if self.max_size < self.min_size:
            raise ValueError('max_size must be greater than min_size')
        # In a simple case, assume we have a function that returns a value between 0.0 and 1.0
        # for any given point along this axis. That value is the normalized finger width that is
        # desired at that location. Ideally, we could go through and build a list of deltas, but
        # we don't know exactly how many are needed, and we can't go below min_size. It is also
        # desirable to return symmetric deltas for symmetric functions.
        #
        # The approach here is to add deltas according to the formula:
        #   delta = min_size + class_scalar[i] * (max_size-min_size)
        # for evenly spaced i until the total length exceeds the required length. Once the number
        # of points is known, create the deltas that fit, and then scale them to fit the size.
        n = 3
        required_length = self.get_total_length()
        while True:
            self.re_seed()
            deltas = [self.get_delta(t) for t in linspace(num=n)]
            current_length = sum(deltas)
            if current_length >= required_length:
                break
            else:
                n += 2
        if current_length > required_length:
            n -= 2
        self.re_seed()
        unscaled_deltas = [self.get_delta(t) for t in linspace(num=n)]
        unscaled_length = sum(unscaled_deltas)
        final_deltas = [unscaled_delta * required_length / unscaled_length for unscaled_delta in unscaled_deltas]

        # Interface with the rest of the world is the member self.locations
        running_total = self.min_inset
        locations = []
        for delta in final_deltas[:-1]: # first and last points are implicit
            running_total += delta
            locations.append(running_total)
        return locations

class ConstantFingerLocations(VariableFingerLocations):
    def get_scalar(self, normalized_location):
        return 0.0

class SquishInFingerLocations(VariableFingerLocations):
    def get_scalar(self, normalized_location):
        return 2.0*abs(normalized_location-0.5)

class SquishOutFingerLocations(SquishInFingerLocations):
    def get_scalar(self, normalized_location):
        return 1.0 - SquishInFingerLocations.get_scalar(self, normalized_location)

class SquishMinFingerLocations(VariableFingerLocations):
    def get_scalar(self, normalized_location):
        return normalized_location

class SquishMaxFingerLocations(SquishMinFingerLocations):
    def get_scalar(self, normalized_location):
        return 1.0 - SquishMinFingerLocations.get_scalar(self, normalized_location)

class SquishMultiFingerLocations(VariableFingerLocations):
    def __init__(self, *args, **kwargs):
        VariableFingerLocations.__init__(self, *args, **kwargs)
        self.children['periods'] = FloatEntry('periods', description = 'Number of periods (repetions of pattern)', default = 2.0, minimum=0, maximum=None)

    def get_scalar(self, normalized_location):
        return math.cos(self.periods*math.pi*normalized_location)**2

class RandomFingerLocations(VariableFingerLocations):
    def __init__(self, *args, **kwargs):
        VariableFingerLocations.__init__(self, *args, **kwargs)
        self.children['periods'] = IntegerEntry('periods', description = 'Random Seed (Enter any number, using the same one will give consistent results)',
                                                default=8341, minimum=0, maximum=100000)

    def re_seed(self):
        random.seed(self.seed)

    def get_scalar(self, normalized_location):
        return random.random()

def finger_factory(finger_type, key, points_from_json = True):
    finger_cls = None
    if finger_type.lower() == 'custom':
        finger_cls = CustomFingerLocations
    if finger_type.lower() == 'squishin':
        finger_cls = SquishInFingerLocations
    if finger_type.lower() == 'squishout':
        finger_cls = SquishOutFingerLocations
    if finger_type.lower() == 'squishmin':
        finger_cls = SquishMinFingerLocations
    if finger_type.lower() == 'squishmax':
        finger_cls = SquishMaxFingerLocations
    if finger_type.lower() == 'squishmulti':
        finger_cls = SquishMultiFingerLocations
    if finger_type.lower() == 'equal':
        finger_cls = ConstantFingerLocations
    if finger_type.lower() == 'random':
        finger_cls = RandomFingerLocations
    if finger_type.lower() == 'new':
        finger_cls = SquishMultiFingerLocations
    if finger_cls is None:
        raise ValueError('Unkown finger type. Received: ' + finger_type)
    return finger_cls(key, points_from_json)

class FingerFactory(EntryFactory):
    def __init__(self, factory_entry, **kwargs):
        EntryFactory.__init__(self, factory_entry = factory_entry,
                              factory_fn = finger_factory,
                              factory_kwargs = kwargs)

#################### END OF FINGERS ##################


class RhinoBoxEntry(ParentEntry):
    def __init__(self, *args, **kwargs):
        ParentEntry.__init__(self, *args, **kwargs)
        for extremum in ALL_EXTREMA:
            self.children[str(extremum)] = FloatEntry(extremum,
                                                      description = '{} point on axis {}'.format(extremum.minmax, extremum.axis),
                                                      default = 0 if extremum.is_min() else 1)

    def draw(self):
        box = RhinoBox(self.value)
        return box.draw()

class BoardEntry(ParentEntry):
    def __init__(self, *args, **kwargs):
        ParentEntry.__init__(self, *args, **kwargs)
        self.children['thickness'] = FloatEntry('thickness', description = 'Thickness of the board',
                                                default = DEFAULT_BOARD_THICKNESS)
        self.children['width'] = FloatEntry('width', description = 'Width of the board',
                                            default = DEFAULT_BOARD_WIDTH)
        self.children['length'] = FloatEntry('length', description = 'Length of the board',
                                             default = DEFAULT_BOARD_WIDTH)
        self.children['normal'] = AxisEntry('normal', description = 'Normal axis to the board (thickness direction)',
                                             default = 'z')
        self.children['grain'] = AxisEntry('grain', description = 'Axis along length of the board (follows grain)',
                                           default = 'x')
        self.children['center'] = PointEntry('center', description = 'Location of center point')
        self.children['is_flipped'] = BoolEntry('is_flipped', description = 'Is board flipped (negative axis normal)', default = False)

    def draw(self):
        box = BasicBoard(**self.value)
        return box.draw()

class CubeEntry(ParentEntry):
    def __init__(self, *args, **kwargs):
        ParentEntry.__init__(self, *args, **kwargs)
        self.children['size'] = PointEntry('size', description = 'Size', default = 150.0)
        self.children['thickness'] = FloatEntry('thickness', description = 'Thickness', default = 4.0)
        self.children['laydown'] = BoolEntry('laydown', description = 'Automatically lay down parts', default = False)
        for extremum in ALL_EXTREMA:
            key = 'enable_{}'.format(extremum)
            self.children[key] = BoolEntry(key,
                                           description = 'Draw {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                           default = True)
        for extremum in ALL_EXTREMA:
            key = 'num_boards_wide_{}'.format(extremum)
            self.children[key] = IntegerEntry(key,
                                              description = 'Number of boards that make up the width of {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                              default = 1)

        for extremum in ALL_EXTREMA:
            key = 'num_boards_long_{}'.format(extremum)
            self.children[key] = IntegerEntry(key,
                                              description = 'Number of boards that make up the length of {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                              default = 1)
        for axis in ALL_AXES:
            key = 'finger_{}'.format(axis)
            type_entry = FingerTypeEntry('type', description = 'axis {}'.format(extremum.axis))
            self.factories[key] = FingerFactory(type_entry, points_from_json = False)

    def is_extremum_enabled(self, extremum):
        return self.children['enable_{}'.format(extremum)].value

    def get_num_boards_wide(self, extremum):
        return self.children['num_boards_wide_{}'.format(extremum)].value

    def get_num_boards_long(self, extremum):
        return self.children['num_boards_long_{}'.format(extremum)].value

    def get_width_axis(self, extremum):
        if extremum.axis == 'x':
            return 'z'
        elif extremum.axis == 'y':
            return 'z'
        elif extremum.axis == 'z':
            return 'y'

    def get_length_axis(self, extremum):
        if extremum.axis == 'x':
            return 'y'
        elif extremum.axis == 'y':
            return 'x'
        elif extremum.axis == 'z':
            return 'x'

    def draw(self):
        thickness = self.thickness
        sizes = self.children['size'].value
        laydown = self.children['laydown'].value
        boards = {}
        board_ids_per_face = {}
        finger_entries = {axis: self.factories['finger_{}'.format(axis)].get_entry() for axis in ALL_AXES}
        curve_ids = {}
        for axis in ALL_AXES:
            min_point = {key: 0.0 for key in ALL_AXES}
            min_point[axis] -= 0.5*sizes[axis]
            max_point = {key: 0.0 for key in ALL_AXES}
            max_point[axis] += 0.5*sizes[axis]
            curve_ids[axis] = rs.AddLine(dict_to_point(min_point), dict_to_point(max_point))
            finger_entries[axis].set_curve(curve_ids[axis])
            finger_entries[axis].intersector_size = 1.01 * max([sizes[other_axis] for other_axis in get_other_axes(axis)])
        for extremum in ALL_EXTREMA:
            if self.is_extremum_enabled(extremum):
                width_axis = self.get_width_axis(extremum)
                length_axis = self.get_length_axis(extremum)
                normal_axis = extremum.axis
                center_point = {'x':0, 'y':0, 'z':0}
                center_point[normal_axis] = 0.5 * extremum.to_minmax_sign() * (sizes[normal_axis] - thickness)
                board = BasicBoard(
                        thickness = thickness,
                        width = sizes[width_axis],
                        length = sizes[length_axis],
                        normal = normal_axis,
                        grain = length_axis,
                        center = center_point,
                        is_flipped = extremum.is_min(),
                        num_boards_wide = self.get_num_boards_wide(extremum),
                        num_boards_long = self.get_num_boards_long(extremum),
                        )
                boards[extremum] = board
                board_ids_per_face[extremum] = board.draw(width_finger_generator = finger_entries[width_axis],
                                                          length_finger_generator = finger_entries[length_axis])

        # Transform to tracking per axis for applying finger joints more efficiently
        board_ids_per_axis = {axis:{} for axis in ALL_AXES}
        for axis in ALL_AXES:
            for minmax in ALL_MINMAX:
                extremum = Extremum(axis, minmax)
                if self.is_extremum_enabled(extremum):
                    board_ids_per_axis[axis][minmax] = board_ids_per_face[extremum]
        for axis in ALL_AXES:
            other_axes = get_other_axes(axis)
            # shallow copies, the dictionaries are modified within the apply function
            primary_ids_dict = board_ids_per_axis[other_axes[1]]
            secondary_ids_dict = board_ids_per_axis[other_axes[0]]
            finger_entries[axis].apply_in_groups(primary_ids_dict, secondary_ids_dict)

        rs.DeleteObjects(curve_ids.values())

        # Transform back to working "per board" if laying down
        for extremum in ALL_EXTREMA:
            if self.is_extremum_enabled(extremum):
                board_ids_per_face[extremum] = board_ids_per_axis[extremum.axis][extremum.minmax]
                if board_ids_per_face[extremum] is None:
                    raise ValueError('Unable to collect object ids of jointed boards')

        # You do not expect laydown to work anymore, because you worked per-axis for fingers instead of per-face
        # One general solution is to extend the finger jointing to groups of groups so that boards can be kept together
        if laydown:
            for extremum in ALL_EXTREMA:
                if self.is_extremum_enabled(extremum):
                    board_ids_per_face[extremum] = boards[extremum].laydown(board_ids_per_face[extremum])
                    translation = [0.0, 0.0, 0.0]
                    if extremum.axis == 'z':
                        translation[1] += extremum.to_minmax_sign() * (0.5*sizes['y'] + 0.5*sizes['z'] + thickness)
                    if extremum.axis == 'y' or extremum.axis == 'z':
                        translation[0] += extremum.to_minmax_sign() * (0.5*sizes['x'] + 1.5*thickness)
                    if extremum.axis == 'x':
                        translation[0] += extremum.to_minmax_sign() * (1.0*sizes['x'] + 0.5*sizes['y'] + 2.5*thickness)
                    board_ids_per_face[extremum] = rs.MoveObjects(board_ids_per_face[extremum], translation)


class TwoPlyCubeEntry(ParentEntry):
    def __init__(self, *args, **kwargs):
        ParentEntry.__init__(self, *args, **kwargs)
        self.num_layers = 2
        self.children['size'] = PointEntry('size', description = 'Size', default = 150.0)
        self.children['thickness'] = FloatEntry('thickness', description = 'Thickness', default = 4.0)
        self.children['laydown'] = BoolEntry('laydown', description = 'Automatically lay down parts', default = False)
        for extremum in ALL_EXTREMA:
            key = 'enable_{}'.format(extremum)
            self.children[key] = BoolEntry(key,
                                           description = 'Draw {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                           default = True)
        for iLayer in range(self.num_layers):
            for extremum in ALL_EXTREMA:
                key = 'num_boards_wide_{}_layer_{}'.format(extremum, iLayer)
                self.children[key] = IntegerEntry(key,
                                                  description = 'Number of boards that make up the width of {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                                  default = 1)

            for extremum in ALL_EXTREMA:
                key = 'num_boards_long_{}_layer_{}'.format(extremum, iLayer)
                self.children[key] = IntegerEntry(key,
                                                  description = 'Number of boards that make up the length of {} face on axis {}'.format(extremum.minmax, extremum.axis),
                                                  default = 1)
        for axis in ALL_AXES:
            key = 'finger_{}'.format(axis)
            type_entry = FingerTypeEntry('type', description = 'axis {}'.format(extremum.axis))
            self.factories[key] = FingerFactory(type_entry, points_from_json = False)

    def is_extremum_enabled(self, extremum):
        return self.children['enable_{}'.format(extremum)].value

    def get_num_boards_wide(self, extremum, iLayer):
        return self.children['num_boards_wide_{}_layer_{}'.format(extremum, iLayer)].value

    def get_num_boards_long(self, extremum, iLayer):
        return self.children['num_boards_long_{}_layer_{}'.format(extremum, iLayer)].value

    def get_sizes(self, iLayer):
        rv = copy.deepcopy(self.size)
        for axis_key in rv:
            rv[axis_key] -= 2 * iLayer * self.thickness
        return rv

    def get_width_axis(self, extremum):
        if extremum.axis == 'x':
            return 'z'
        elif extremum.axis == 'y':
            return 'z'
        elif extremum.axis == 'z':
            return 'y'

    def get_length_axis(self, extremum):
        if extremum.axis == 'x':
            return 'y'
        elif extremum.axis == 'y':
            return 'x'
        elif extremum.axis == 'z':
            return 'x'

    def get_board_key(self, extremum, iLayer):
        return '{}_{}'.format(extremum, iLayer)

    def get_axis_key(self, axis, iLayer):
        return '{}_{}'.format(axis, iLayer)

    def draw(self):
        thickness = self.thickness
        laydown = self.children['laydown'].value
        boards = {}
        board_ids_per_face = {}
        finger_entries = {axis: self.factories['finger_{}'.format(axis)].get_entry() for axis in ALL_AXES}
        curve_ids = {}
        for axis in ALL_AXES:
            min_point = {key: 0.0 for key in ALL_AXES}
            min_point[axis] -= 0.5*self.size[axis]
            max_point = {key: 0.0 for key in ALL_AXES}
            max_point[axis] += 0.5*self.size[axis]
            curve_ids[axis] = rs.AddLine(dict_to_point(min_point), dict_to_point(max_point))
            finger_entries[axis].set_curve(curve_ids[axis])
            finger_entries[axis].intersector_size = 1.01 * max([self.size[other_axis] for other_axis in get_other_axes(axis)])

        for iLayer in range(self.num_layers):
            sizes = self.get_sizes(iLayer)
            for extremum in ALL_EXTREMA:
                if self.is_extremum_enabled(extremum):
                    width_axis = self.get_width_axis(extremum)
                    length_axis = self.get_length_axis(extremum)
                    normal_axis = extremum.axis
                    center_point = {'x':0, 'y':0, 'z':0}
                    center_point[normal_axis] = 0.5 * extremum.to_minmax_sign() * (sizes[normal_axis] - thickness)
                    board = BasicBoard(
                            thickness = thickness,
                            width = self.size[width_axis],
                            length = self.size[length_axis],
                            normal = normal_axis,
                            grain = length_axis,
                            center = center_point,
                            is_flipped = extremum.is_min(),
                            num_boards_wide = self.get_num_boards_wide(extremum, iLayer),
                            num_boards_long = self.get_num_boards_long(extremum, iLayer),
                            )
                    board_key = self.get_board_key(extremum, iLayer)
                    boards[board_key] = board
                    board_ids_per_face[board_key] = board.draw(width_finger_generator = finger_entries[width_axis],
                                                               length_finger_generator = finger_entries[length_axis])

        board_ids_per_axis = {}
        # Transform to tracking per axis for applying finger joints more efficiently
        for iLayer in range(self.num_layers):
            for axis in ALL_AXES:
                board_ids_per_axis[self.get_axis_key(axis, iLayer)] = {}
        for iLayer in range(self.num_layers):
            for axis in ALL_AXES:
                for minmax in ALL_MINMAX:
                    extremum = Extremum(axis, minmax)
                    if self.is_extremum_enabled(extremum):
                        axis_key = self.get_axis_key(axis, iLayer)
                        board_key = self.get_board_key(extremum, iLayer)
                        board_ids_per_axis[axis_key][minmax] = board_ids_per_face[board_key]
        # First apply fingers for the inner layer
        for iLayer in range(self.num_layers):
            for axis in ALL_AXES:
                other_axes = get_other_axes(axis)
                # shallow copies, the dictionaries are modified within the apply function
                primary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[1], iLayer)]
                secondary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[0], iLayer)]
                finger_entries[axis].apply_in_groups(primary_ids_dict, secondary_ids_dict, polarity = iLayer%2)
        # Second push every third inner finger through the outer layer
        for iLayer in [0, ]:
            for axis in ALL_AXES:
                other_axes = get_other_axes(axis)
                # shallow copies, the dictionaries are modified within the apply function
                primary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[1], iLayer+1)]
                secondary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[0], iLayer)]
                finger_entries[axis].apply_in_groups(secondary_ids_dict, primary_ids_dict, draw_type = 'mini', polarity=0)
            for axis in ALL_AXES:
                other_axes = get_other_axes(axis)
                # shallow copies, the dictionaries are modified within the apply function
                primary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[1], iLayer)]
                secondary_ids_dict = board_ids_per_axis[self.get_axis_key(other_axes[0], iLayer+1)]
                finger_entries[axis].apply_in_groups(primary_ids_dict, secondary_ids_dict, draw_type = 'mini_offset')

        rs.DeleteObjects(curve_ids.values())

        # Transform back to working "per board" if laying down
        for iLayer in range(self.num_layers):
            # Transform back to working "per board" if laying down
            for extremum in ALL_EXTREMA:
                if self.is_extremum_enabled(extremum):
                    axis_key = self.get_axis_key(extremum.axis, iLayer)
                    board_key = self.get_board_key(extremum, iLayer)
                    board_ids_per_face[board_key] = board_ids_per_axis[axis_key][extremum.minmax]


        if laydown:
            for iLayer in range(self.num_layers):
                sizes = self.size
                for extremum in ALL_EXTREMA:
                    if self.is_extremum_enabled(extremum):
                        board_key = self.get_board_key(extremum, iLayer)
                        board_ids_per_face[board_key] = boards[board_key].laydown(board_ids_per_face[board_key])
                        translation = [0.0, 0.0, 0.0]
                        if extremum.axis == 'z':
                            translation[1] += extremum.to_minmax_sign() * (0.5*sizes['y'] + 0.5*sizes['z'] + thickness)
                        if extremum.axis == 'y' or extremum.axis == 'z':
                            translation[0] += extremum.to_minmax_sign() * (0.5*sizes['x'] + 1.5*thickness)
                        if extremum.axis == 'x':
                            translation[0] += extremum.to_minmax_sign() * (1.0*sizes['x'] + 0.5*sizes['y'] + 2.5*thickness)
                        translation[1] += (sizes['z'] + 2.0*sizes['y'] + thickness) * (self.num_layers - iLayer - 1)
                        print board_ids_per_face[board_key]
                        board_ids_per_face[board_key] = rs.MoveObjects(board_ids_per_face[board_key], translation)


def test_reload():
    entries = []
    entries.append(FingerTypeEntry('finger_type'))
    entries.append(FloatEntry(
        key = 'min_joint_size',
        description = 'Minimum joint size',
        default = 2,
        minimum = 0))
    entries.append(RhinoBoxEntry('something'))

    copied_entries = [entry for entry in entries]

    for entry in entries:
        entry.from_user_prompt()

    config = {}
    for entry in entries:
        entry.to_property_tree(config)

    config_str = json.dumps(config, sort_keys=True)
    #config_str = json.dumps(config, sort_keys=True, indent=4)
    print(config_str)

    reloaded_config = json.loads(config_str)
    for entry in copied_entries:
        entry.from_property_tree(reloaded_config)
        print(entry)

def test_draw_something():
    entry = RhinoBoxEntry('example_board')
    entry_json = '{"example_board": {"x_max": 1000.0, "x_min": 0.0, "y_max": 500.0, "y_min": 0.0, "z_max": 700.0, "z_min": 0.0}}'
    loaded_config = json.loads(entry_json)
    entry.from_property_tree(loaded_config)
    entry.draw()

def test_draw_board():
    entry = BoardEntry('example_board')
    if 1:
        entry_json = '{"example_board": {"thickness": 4.0, "width": 150.0, "length": 300.0, "normal": "y", "grain": "z", "center": {"x": -200.0, "y": 100.0, "z": 20.0}, "is_flipped": false}}'
        loaded_config = json.loads(entry_json)
        entry.from_property_tree(loaded_config)
    else:
        entry.from_user_prompt()
    config = {}
    entry.to_property_tree(config)
    config_str = json.dumps(config)
    print(config_str)
    entry.draw()

def test_draw_cube():
    entry = CubeEntry('example_cube')
    if 1:
        entry_json = '''
{
    "example_cube": {
        "size": {
            "x": 250.0,
            "y": 150.0,
            "z": 75.0
        },
        "thickness": 4.0,
        "laydown": true,
        "enable_x_min": true,
        "enable_x_max": true,
        "enable_y_min": true,
        "enable_y_max": true,
        "enable_z_min": false,
        "enable_z_max": false,
        "num_boards_wide_x_min": 1,
        "num_boards_wide_x_max": 1,
        "num_boards_wide_y_min": 2,
        "num_boards_wide_y_max": 1,
        "num_boards_wide_z_min": 2,
        "num_boards_wide_z_max": 1,
        "num_boards_long_x_min": 1,
        "num_boards_long_x_max": 2,
        "num_boards_long_y_min": 2,
        "num_boards_long_y_max": 1,
        "num_boards_long_z_min": 1,
        "num_boards_long_z_max": 1,
        "finger_x": {
            "min_size": 2.5,
            "min_inset": 4.0,
            "max_inset": 4.0,
            "max_size": 6.5,
            "periods": 2.0,
            "type": "SquishMulti"
        },
        "finger_y": {
            "min_size": 2.5,
            "min_inset": 4.0,
            "max_inset": 4.0,
            "max_size": 6.5,
            "periods": 1.0,
            "type": "SquishMulti"
        },
        "finger_z": {
            "min_size": 2.5,
            "min_inset": 4.0,
            "max_inset": 4.0,
            "max_size": 6.5,
            "periods": 1.0,
            "type": "SquishMulti"
        }
    }
}
'''
        loaded_config = json.loads(entry_json)
        entry.from_property_tree(loaded_config)
    else:
        entry.from_user_prompt()
    config = {}
    entry.to_property_tree(config)
    config_str = json.dumps(config, indent=4)
    print(config_str)
    entry.draw()

    # Don't control-z past this point
def test_draw_twoply_cube():
    entry = TwoPlyCubeEntry('example_cube')
    if 1:
        entry_json = '''
{
    "example_cube": {
        "size": {
            "x": 336.5,
            "y": 208.0,
            "z": 180.0
        },
        "thickness": 4.0,
        "laydown": true,
        "enable_x_min": true,
        "enable_x_max": true,
        "enable_y_min": true,
        "enable_y_max": true,
        "enable_z_min": true,
        "enable_z_max": false,
        "num_boards_wide_x_min_layer_0": 1,
        "num_boards_wide_x_max_layer_0": 1,
        "num_boards_wide_y_min_layer_0": 1,
        "num_boards_wide_y_max_layer_0": 1,
        "num_boards_wide_z_min_layer_0": 1,
        "num_boards_wide_z_max_layer_0": 1,
        "num_boards_long_x_min_layer_0": 1,
        "num_boards_long_x_max_layer_0": 1,
        "num_boards_long_y_min_layer_0": 1,
        "num_boards_long_y_max_layer_0": 1,
        "num_boards_long_z_min_layer_0": 1,
        "num_boards_long_z_max_layer_0": 1,
        "num_boards_wide_x_min_layer_1": 1,
        "num_boards_wide_x_max_layer_1": 1,
        "num_boards_wide_y_min_layer_1": 1,
        "num_boards_wide_y_max_layer_1": 1,
        "num_boards_wide_z_min_layer_1": 1,
        "num_boards_wide_z_max_layer_1": 1,
        "num_boards_long_x_min_layer_1": 1,
        "num_boards_long_x_max_layer_1": 1,
        "num_boards_long_y_min_layer_1": 1,
        "num_boards_long_y_max_layer_1": 1,
        "num_boards_long_z_min_layer_1": 1,
        "num_boards_long_z_max_layer_1": 1,
        "finger_x": {
            "min_size": 5.5,
            "min_inset": 8.0,
            "max_inset": 8.0,
            "max_size": 10.5,
            "periods": 2.0,
            "type": "SquishMulti"
        },
        "finger_y": {
            "min_size": 5.5,
            "min_inset": 8.0,
            "max_inset": 8.0,
            "max_size": 10.5,
            "periods": 1.0,
            "type": "SquishMulti"
        },
        "finger_z": {
            "min_size": 5.5,
            "min_inset": 8.0,
            "max_inset": 8.0,
            "max_size": 10.5,
            "periods": 1.0,
            "type": "SquishMulti"
        }
    }
}
'''
        loaded_config = json.loads(entry_json)
        entry.from_property_tree(loaded_config)
    else:
        entry.from_user_prompt()
    config = {}
    entry.to_property_tree(config)
    config_str = json.dumps(config, indent=4)
    print(config_str)
    entry.draw()

def test_apply_fingers():
    entry = finger_factory('SquishMulti', 'finger_style', points_from_json=False)
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
        loaded_config = json.loads(entry_json)
        entry.from_property_tree(loaded_config)
    else:
        entry.from_user_prompt()
    config = {}
    entry.to_property_tree(config)
    config_str = json.dumps(config, indent=4)
    print(config_str)
    entry.init_with_user_curve()
    entry.apply_user()


def main():
    #test_reload()
    #test_draw_something()
    #test_draw_board()
    #test_draw_cube()
    test_draw_twoply_cube()
    #test_apply_fingers()

if __name__ == '__main__':
    main()
