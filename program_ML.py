# _*_ coding: utf-8 _*_
'''
@author ML
reconstruct my project
2018,5,31
'''
from hp import *
import xml
import xml.etree.ElementTree as et
from sympy.geometry import *
import matplotlib.pyplot as plt
from sympy import pi


class dir_seg(object):
    def __init__(self, p1, p2):
        self.p1, self.p2 = p1, p2
        if p1.x == p2.x:
            self.horizontal = False
            self.vertical = True
        elif p1.y == p2.y:
            self.horizontal = True
            self.vertical = False
        self.__init()
    def __init(self):
        p1, p2 = self.p1, self.p2
        self.line = Line(p1, p2)
        self.seg = Segment(p1, p2)
        self.ray = Ray(p1, p2)
        r = Ray(p1, p2)
        p = Point(r.p2 - r.p1)
        len = self.seg.length
        self.normal = Ray((0, 0), (p.y / len, -p.x / len))
        self.dir = Ray((0, 0), (p.x / len, p.y / len))
    def cross_dir(self, seg):
        a = ((self.p2.x-self.p1.x),(self.p2.y - self.p1.y))
        b = ((seg.p2.x-seg.p1.x),(seg.p2.y - seg.p1.y))
        direc = (a[0] * b[1] - b[0] * a[1])
        if direc != 0:
            return direc/abs(direc)
        else:
            return 0
    def get_p1_p2_from_normal(normal, p1, p2):
        r = dir_seg(p1, p2)
        if r.normal.angle_between(normal) == 0:
            return p1, p2
        elif r.normal.angle_between(normal) == pi:
            return p2, p1
        else:
            assert True, "法线设置有误"
    
class ml_boundary():
    def __init__(self, *args, **kwargs):
        polygon_tmp = Polygon(*args, **kwargs)
        self.polygon = Polygon(*polygon_tmp.vertices)
        self.seg_list = []
        self.init_seg()
    def init_seg(self):
        self.seg_list.clear()
        v = self.polygon.vertices
        for i in range(0,len(v) - 1):
            self.seg_list.append(dir_seg(v[i], v[i + 1]))
        self.seg_list.append(dir_seg(v[len(v) - 1],v[0]))
    def draw_b(self,*args):
        if args==():
            type_line = '-'
            C = '#000000'
        else:
            type_line = '--'
            C = '#990033'
        for s in self.seg_list:
            x=[s.p1[0] , s.p2[0]]
            y=[s.p1[1] , s.p2[1]]
            plt.plot(x,y,type_line,color= C )
def xml_read(file):
    boundary_list = []
    xml_file = et.parse(file)
    root =  xml_file.getroot()
    house = xml_file.findall('House')[0]
    fp = house.findall('FloorPlan')[0]
    son = fp.getchildren()
    boundary = son[0].attrib['boundary']
    boundary_list.append(boundary)
    return boundary_list[0]
def set_boundary(bd):
    tw = bd.split(';')
    p_list = []
    for t in tw:
        numl = t[1:-1].split(',')
        p_list.append(Point2D(int(numl[0]),int(numl[1])))
    bdy = ml_boundary(*p_list)
    return bdy

def main(file):
    boundary_date = xml_read(file)
    boundary = set_boundary(boundary_date)
    boundary.draw_b()

    vb_list = get_n_virtual_boundary(boundary)
    for vb in vb_list:
        vb.draw_b(1)
    # vb = get_virtual_boundary(boundary)
    # vb.draw_b(1)
    plt.show()

if __name__ == '__main__':

    file = r'C:\Users\dyrs-ai-win10\Desktop\xmltest\fenquao.xml'
    main(file)