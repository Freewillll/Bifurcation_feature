import sys
from swc_handler import *
from file_io import *
from morph_topo.morphology_features import *
from path_util import *


class Bifurcation(object):

    def __init__(self, idx, morph, line_length=8.0, z_factor=1):

        self.morph = morph
        self.z_factor = z_factor
        topo_tree, seg_dict = self.morph.convert_to_topology_tree()
        self.topo = Topology(topo_tree)
        self.seg_dict = seg_dict
        self.line_length = line_length
        self.bif_idx = idx
        self.level = self.topo.order_dict[self.bif_idx]

    def get_angles(self):
        child_1_idx = self.morph.child_dict[self.bif_idx][0]
        child_2_idx = self.morph.child_dict[self.bif_idx][1]
        par_idx = self.morph.pos_dict[self.bif_idx][6]

        c = np.array(self.morph.pos_dict[self.bif_idx][2:5])

        c1 = find_point_by_distance(pt=c, anchor_idx=child_1_idx, is_parent=False, morph=self.morph,
                                    dist=self.line_length, return_center_point=True, spacing=(1,1,self.z_factor))
        c2 = find_point_by_distance(pt=c, anchor_idx=child_2_idx, is_parent=False, morph=self.morph,
                                    dist=self.line_length, return_center_point=True, spacing=(1,1,self.z_factor))
        p1 = find_point_by_distance(pt=c, anchor_idx=par_idx, is_parent=True, morph=self.morph,
                                    dist=self.line_length, return_center_point=True, spacing=(1,1,self.z_factor))

        ang_c1_p1 = calc_included_angles_from_coords(c, c1, p1, return_rad=True)[0]
        ang_c2_p1 = calc_included_angles_from_coords(c, c2, p1, return_rad=True)[0]
        ang_c1_c2 = calc_included_angles_from_coords(c, c1, c2, return_rad=True)[0]

        ang_c1_p1, ang_c2_p1 = sorted([ang_c1_p1, ang_c2_p1])
        return ang_c1_p1, ang_c2_p1, ang_c1_c2


def bif_feature(tree, line_length):
    feature = []
    morph = Morphology(tree)
    morph.get_critical_points()
    print(len(morph.bifurcation))
    for idx in morph.bifurcation:
        bif = Bifurcation(idx, morph, line_length)
        ang_c1_p1, ang_c2_p1, ang_c1_c2 = bif.get_angles()
        print(f'ang1:{ang_c1_p1}, ang2:{ang_c2_p1}, ang3:{ang_c1_c2}, level:{bif.level}')
        feature.append((ang_c1_p1, ang_c2_p1, ang_c1_c2, bif.level))
    return feature


if __name__ == '__main__':
    swcfile = 'C:/Users/fmind/Desktop/17109_1701_x8048_y22277.semi_r.swc'
    tree = parse_swc(swcfile)
    fea = bif_feature(tree, 10)
