# _*_ coding: utf-8 _*_
from program_ML import *
import math


def get_adj_seg(ver, boundary):
    '''返回同一点的相邻边'''
    tp_adj_seg_list = []
    for seg in boundary.seg_list:
        if seg.seg.contains(ver):
            tp_adj_seg_list.append(seg)
    return tp_adj_seg_list


def another_p(line, point):
    '''返回一条边的另一个顶点'''
    if line.p1 == point:
        return line.p2
    else:
        return line.p1


def point_distance(p1, p2):
    a = (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2
    a = math.sqrt(a)
    return a


def get_inner_line_list(boundary, concave_list=[]):
    ''' 得到向内凹的线'''
    if concave_list == []:
        concave_list, _ = concave_or_convex(boundary)
    inner_line_l = []
    for i in boundary.seg_list:
        if i.p1 in concave_list and i.p2 in concave_list:
            inner_line_l.append(i)
    return inner_line_l


def concave_or_convex_point(la, lb):
    return la.cross_dir(lb)


def concave_or_convex(boundary):
    '''concave表示凹点，convex表示凸点
    返回两个列表
    '''
    concave_list = []
    convex_list = []
    for v in boundary.polygon.vertices:
        v_list = get_adj_seg(v, boundary)
        a = another_p(v_list[0], v)
        if a == v_list[0].p1:
            la = v_list[0]
            lb = v_list[1]
        else:
            la = v_list[1]
            lb = v_list[0]
        tag = concave_or_convex_point(la, lb)
        if tag == 1:
            # 凹点
            concave_list.append(v)
        else:  # tag==-1:
            # 凸出点
            convex_list.append(v)
    return concave_list, convex_list


def cover_range(seg):
    if seg.p1.x == seg.p2.x:
        if seg.p1.y > seg.p2.y:
            return (seg.p2.y, seg.p1.y)
        else:
            return (seg.p1.y, seg.p2.y)
    else:
        if seg.p1.x > seg.p2.x:
            return (seg.p2.x, seg.p1.x)
        else:
            return (seg.p1.x, seg.p2.x)


def ray_ext(a, seg_list):
    vecter = a.dir.p2
    a2 = a.p2
    ray_line = dir_seg(a2, a2 + vecter * 100000000)
    ext_list = []
    for i in seg_list:
        if ray_line.seg.intersection(i.seg) != []:
            ext_list.append((i.line.distance(a2), ray_line.seg.intersection(i.seg)[0], i))
            # 返回元组有三项，距离，交点，交点线
    if ext_list == []:
        raise Exception('传入射线和列表有问题', a.p1, a.p2, vecter)
    else:
        ext_list.sort(key=lambda x: x[0])
    if ext_list[0][0] != 0:
        raise Exception('射线起始点不在凹凸点上，这个图形有问题', a.p1, a.p2)
    return ext_list


def is_staggerd(seg, area_range, need_tag=0):
    '''  测试边， 测试范围,测试范围自己先提前调用cover——range得到
    暂时返回值只用到0,1，但是用标记表明了各种相交情况，当需要时就可以更改返回值使用
    需要标记时第三个参数设置1，返回标志
    '''
    p_l = cover_range(seg)
    ar = [area_range[0], area_range[1]]
    ar.sort()
    if p_l[1] <= ar[0] or p_l[0] >= ar[1]:
        # waimian
        tag = 'n'
        re = 0
    elif p_l[0] <= ar[0] and p_l[1] > ar[0] and p_l[1] < ar[1]:
        tag = 'l'
        re = 1
    elif p_l[0] > ar[0] and p_l[1] < ar[1]:
        tag = 'm'
        re = 1
    elif p_l[0] > ar[0] and p_l[0] < ar[1] and p_l[1] > ar[1]:
        tag = 'r'
        re = 1
    else:
        tag = 'a'  # all
        re = 1
    if need_tag:
        return tag
    else:
        return re


def squar(a, b):
    '''对角点，算面积'''
    return abs((a.x - b.x) * (a.y - b.y))


def avoid_func(seg, avo_list):
    ''' 一条边与列表中所有边，在某个坐标方向上，坐标区域的重合关系'''
    n, a, l, r, m = 0, 0, 0, 0, 0
    s_range = cover_range(seg)
    lp, rp = None, None
    m_range = []

    for i in avo_list:
        ran = cover_range(i)
        if ran[0] <= s_range[0] and s_range[0] < ran[1] and ran[1] < s_range[1]:
            l += 1
            lp = ran[1]
        elif ran[0] > s_range[0] and ran[0] < s_range[1] and ran[1] >= s_range[1]:
            r += 1
            rp = ran[0]
        elif s_range[0] < ran[0] and s_range[1] > ran[1]:
            m += 1
            m_range.append(ran)
        elif s_range[0] >= ran[0] and s_range[1] <= ran[1]:
            a += 1
        else:
            n += 1
    return [n, a, l, r, m], [lp, rp, m_range]


def cut_line(a, b, l):
    ''' 一条线，沿线方向上两个位置顶点切割'''
    if l.p1.y == l.p2.y:
        # 水平线
        y = l.p1.y
        return dir_seg(Point2D(a, y), Point2D(b, y))

    else:
        x = l.p1.x
        return dir_seg(Point2D(x, a), Point2D(x, b))


def re_shape(a, c):
    '''对角线两个点，重构四边形'''
    if a.x < c.x:
        a, c = c, a
    if a.y > c.y:
        b = Point2D(a.x, c.y)
        d = Point2D(c.x, a.y)
    else:
        b = Point2D(c.x, a.y)
        d = Point2D(a.x, c.y)

    return [a, b, c, d]


def line_ext_max_final(l, praline, boundary, concave_list):
    inside_list = []
    a_dict = {}
    global vector, np_poi
    '''同侧化praline里面的规避边'''
    # 水平线
    vector = l.normal.p2[1]
    np_poi = l.p1.y
    x_rang = cover_range(l)
    for pl in praline:
        ai = pl.p1.y - l.p1.y
        aai = abs(ai)
        if ai / aai == l.normal.p2[1] and is_staggerd(pl, x_rang):
            if a_dict.get(aai) == None:
                a_dict[aai] = []
            a_dict[aai].append(pl)

    inside_list = [(x1, a_dict[x1]) for x1 in a_dict.keys()]
    inside_list.sort(key=lambda x: x[0])

    # d_list全局答案存放
    global d_list
    d_list = []

    def max_rec(tl, inside_list):
        in_list = inside_list.copy()
        try:
            il = in_list[0]
        except:
            raise Exception('矩形曼延处，曼延线没有平行且坐标区域重合的线，区域不闭合', tl.p1, tl.p2)
        tl_range = cover_range(tl)
        true_p = np_poi + vector * il[0]
        in_list.__delitem__(0)
        af, ap = avoid_func(tl, il[1])
        if af[1] or af[2] or af[3] or af[4]:
            s = il[0] * tl.seg.length
            d_list.append((s, cover_range(tl), np_poi, true_p))
        else:
            max_rec(tl, in_list)
        '''更新下波数据,tl要更新，如果要进入递归，则要转换相应参数并且结束这个for循环'''
        if not af[1]:
            # m的情况
            cp = [tl_range[0], tl_range[1]]
            m_range = ap[2]
            for mr in m_range:
                cp.append(mr[0])
                cp.append(mr[1])
            if af[2]:
                cp.append(ap[0])
                cp.sort()
                cp.__delitem__(0)
            if af[3]:
                cp.append(ap[1])
                cp.sort()
                cp.__delitem__(-1)
            cp.sort()
            l = len(cp)
            if l % 2 == 1:
                raise Exception('位于矩形曼延处代码，检查是否有凹陷的单点线，\
                                凹陷处应该为墙而不是单线段', a.p1, a.p2)

            l = int(l / 2)
            for i in range(l):
                tl = cut_line(cp[2 * i], cp[2 * i + 1], tl)
                max_rec(tl, in_list)

    max_rec(l, inside_list)

    d_list.sort(key=lambda x: x[0], reverse=True)
    max_s = d_list[0][0]
    x_rang = d_list[0][1]
    y_rang = (d_list[0][2], d_list[0][3])
    a = Point2D(x_rang[0], y_rang[0])
    c = Point2D(x_rang[1], y_rang[1])

    p_list = [a, c]

    return max_s, p_list


def largest_trangle_final(vex, ho_line, boundary, concave_list=[]):
    horizontal_line_list = [x for x in boundary.seg_list if x.horizontal == True]
    vertical_list = [x for x in boundary.seg_list if x.vertical == True]
    an_list = []
    f_ip = []
    ip = [ho_line.p1, ho_line.p2]
    for p in range(2):
        a, b = ip[p], ip[1 - p]
        ab = dir_seg(a, b)
        while b in concave_list:
            b_l = ray_ext(ab, vertical_list)
            b = b_l[1][1]
            t_dir_l = dir_seg(b, another_p(b_l[1][2], b))
            ab = dir_seg(a, b)

            if t_dir_l.dir.p2 == ho_line.normal.p2:
                break

        f_ip.append(b)
    a, b = f_ip[0], f_ip[1]
    p1, p2 = dir_seg.get_p1_p2_from_normal(ho_line.normal, a, b)
    l = dir_seg(p1, p2)
    s, p_list = line_ext_max_final(l, horizontal_line_list, boundary, concave_list)
    an_list.append((s, p_list))
    an_list.sort(key=lambda x: x[0], reverse=True)
    return an_list[0][0], an_list[0][1]


def get_virtual_boundary(boundary):
    '''新版本的虚拟边界提取函数'''
    concave_list, convex_list = concave_or_convex(boundary)
    a_list = []

    # x_vec_line = dir_seg(Point2D(0, 0), Point2D(1, 0))
    # horizontal_line_list = get_paralleled_line(x_vec_line, boundary, dir_seg)
    horizontal_line_list = [x for x in boundary.seg_list if x.horizontal == True]
    horizontal_line = [x for x in boundary.seg_list if (x.horizontal == True and x.normal.p2.y > 0)]

    for se in horizontal_line:
        vex = se.p1
        s, p_list = largest_trangle_final(vex, se, boundary, concave_list)
        a_list.append((s, p_list))
    a_list.sort(key=lambda x: x[0], reverse=True)
    p_list = re_shape(a_list[0][1][0], a_list[0][1][1])
    virtual_boundary = ml_boundary(*[p_list[0], p_list[1], p_list[2], p_list[3]])

    return virtual_boundary


def get_virtual_boundary_all_rec(boundary):
    '''新版本的虚拟边界提取函数'''
    concave_list, convex_list = concave_or_convex(boundary)
    a_list = []

    # x_vec_line = dir_seg(Point2D(0, 0), Point2D(1, 0))
    # horizontal_line_list = get_paralleled_line(x_vec_line, boundary, dir_seg)
    horizontal_line = [x for x in boundary.seg_list if x.horizontal == True]

    for se in horizontal_line:
        vex = se.p1
        s, p_list = largest_trangle_final(vex, se, boundary, concave_list)
        a_list.append((s, p_list))
    a_list.sort(key=lambda x: x[0], reverse=True)

    return a_list[0]


def get_vector_seg(ver, boundary):
    for i in boundary.seg_list:
        if i.p1 == ver:
            return i


def get_point_belong_seg(ver, boundary):
    seg_list = []
    for seg in boundary.seg_list:
        if seg.seg.contains(ver):
            seg_list.append(seg)
    return seg_list


def divide_rec(adj_pl, not_adjpl, boundary, virtual_boundary, v_b_list=[]):
    '''考虑到第二重深度，直接用第一层的vb, 暂不考虑多个v_b_list
    因为墙是有厚度的，暂时不认为会有单条向内凹的线，即不会有要分的区是以一条边隔开的
    '''

    def seg_contain(seg, poi_list):
        con_list = []
        for poi in poi_list:
            if seg.seg.contains(poi):
                con_list.append(poi)
        return con_list

    vb_point_l = virtual_boundary.polygon.vertices
    adj_pl.extend([x for x in vb_point_l if x not in adj_pl])
    roi_list = []
    while not_adjpl != []:
        sp = not_adjpl[0]
        not_adjpl.remove(sp)
        ver0, ver1 = sp, None
        vec_line = get_vector_seg(ver0, boundary)
        ver1 = vec_line.p2
        roi_p_list = [ver0]
        while ver1 != sp:

            '''应该是先得到所在的向量线，该线沿边界走到虚拟边区，交点会有几个情况，
                    1，没有交点，最正常的走法
                    2，垂直相交，直接走到，向量线p2赋给v1直接走到v1并考虑进入虚边走法
                    3，共线相交，
                        3.1考虑p2点跑到虚拟边界上,可能在虚边界刚好的某个顶点上
                        3.2 p2跑超过了虚拟边界
                    不会是p1点交虚边界，p1交就会是上个循环的p2交，p1,p2同时在虚拟边界上的跑法是错误的
            
                射线在虚边界延伸需要另一套走法
                    1，走到虚边界点，
                    2，走的虚边界垂直实际边界
                    3，共线
            '''

            judge_list = seg_contain(vec_line, adj_pl)
            if judge_list == []:
                # 更新下一轮
                ver0 = ver1
                vec_line = get_vector_seg(ver0, boundary)
                ver1 = vec_line.p2

            elif ver0 not in adj_pl:

                if len(judge_list) == 1:
                    # 所在虚边界的p1和该点重构一条线段
                    ver0 = judge_list[0]
                    t_vector_line_list = get_point_belong_seg(ver0, virtual_boundary)
                    if len(t_vector_line_list) == 1:
                        t_vector_line = t_vector_line_list[0]
                    else:
                        # 刚好在角点
                        for seg in t_vector_line_list:
                            if seg.cross_dir(vec_line) == 0:
                                continue
                            else:
                                t_vector_line = seg
                                break
                    ver1 = t_vector_line.p1
                    vec_line = dir_seg(ver0, ver1)

                elif len(judge_list) >= 2:
                    judge_list.sort(key=lambda x: point_distance(ver0, x))
                    ver0 = judge_list[0]
                    t = ver1
                    t_vector_line_list = get_adj_seg(ver0, virtual_boundary)
                    for seg in t_vector_line_list:
                        if seg.p2 == ver0:
                            ver1 = seg.p1
                    if t == ver1:
                        raise Exception('代码让图形走错路了')
                    vec_line = dir_seg(ver0, ver1)

                else:
                    raise Exception('理论不可能，线不会从虚拟区域中穿过去的，检查一级分区')
            elif ver0 in adj_pl:
                if len(judge_list) == 1:
                    ver0 = ver1
                    vec_line = get_vector_seg(ver0, boundary)
                    ver1 = vec_line.p2
                else:

                    judge_list.sort(key=lambda x: point_distance(ver0, x))
                    # 第一点肯定是ver0本身
                    ver1 = judge_list[1]
                    t_l = dir_seg(ver0, ver1)
                    '''t_p是出去射线源点外最近的点，他可能有很多情况，将他重构成一条射线，判断延伸点情况'''
                    t_vector_line_list = get_adj_seg(ver1, boundary)
                    if t_vector_line_list == []:
                        v_line_list = get_adj_seg(ver1, virtual_boundary)
                        for seg in v_line_list:
                            if seg.line.contains(ver0):
                                continue
                            else:
                                t_vector_line = seg
                                break
                        ver0 = ver1
                        if t_vector_line.p1 == ver1:
                            ver1 = t_vector_line.p2
                        else:
                            ver1 = t_vector_line.p1
                        vec_line = dir_seg(ver0, ver1)
                    elif len(t_vector_line_list) == 1:
                        ver0 = ver1
                        ver1 = t_vector_line_list[0].p2
                        vec_line = dir_seg(ver0, ver1)

                    elif len(t_vector_line_list) == 2:
                        for seg in t_vector_line_list:
                            if seg.cross_dir(t_l) == 0:
                                continue
                            else:
                                vec_line = seg
                                break
                        ver0 = ver1
                        ver1 = vec_line.p2
                        vec_line = dir_seg(ver0, ver1)
                    else:
                        raise Exception('一点不该对应超过两条线，线段错误')
            else:
                raise Exception('unknow')
            try:
                if ver0 in not_adjpl:
                    not_adjpl.remove(ver0)
            except:
                raise Exception('矩形分割处，出现该异常应为区域边界是否有隐含点未被删除')
            roi_p_list.append(ver0)
        bd = ml_boundary(*roi_p_list)
        roi_list.append(bd)
    return roi_list


def get_n_virtual_boundary(boundary, virtual_boundary_list=[], n=2):
    '''调用这个函数就会递归的调用函数本身，直到得到第n级别的提取矩形
        在设计函数实现时暂时未考虑第三重即以上递归调用时虚边界的影响，固现在默认n=2
    '''

    n_count = len(virtual_boundary_list)
    if n_count == 0:
        virtual_boundary = get_virtual_boundary(boundary)
        virtual_boundary_list.append(virtual_boundary)
    else:
        virtual_boundary = virtual_boundary_list[-1]
    adj_point_list = []
    not_adj_point_list = []
    for seg in virtual_boundary.seg_list:
        for i in boundary.polygon.vertices:
            if seg.seg.contains(i) and i not in adj_point_list:
                adj_point_list.append(i)
    for i in boundary.polygon.vertices:
        if i not in adj_point_list:
            not_adj_point_list.append(i)

    # 调用分区函数
    roi_list = divide_rec(adj_point_list, not_adj_point_list, boundary, virtual_boundary, virtual_boundary_list)

    a_list = []
    for roi in roi_list:
        v_list = get_virtual_boundary_all_rec(roi)
        a_list.append(v_list)
    a_list.sort(key=lambda x: x[0], reverse=True)

    p_list = re_shape(a_list[n_count][1][0], a_list[n_count][1][1])

    virtual_boundary = ml_boundary(*[p_list[0], p_list[1], p_list[2], p_list[3]])
    virtual_boundary_list.append(virtual_boundary)

    n_count = len(virtual_boundary_list)
    if n_count >= n:
        # 退出递归-
        return virtual_boundary_list

        # virtual_boundary_list = get_n_virtual_boundary(boundary, n, virtual_boundary_list)
        # return virtual_boundary_list
