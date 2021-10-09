import bpy
import gpu
import bgl
from gpu_extras.batch import batch_for_shader
from mathutils import Vector
from mathutils import Quaternion
from dataclasses import dataclass
from typing import List
import time

# Vector @ Vector = dot
# Quaternion @ Vector = rotate

# s = obj.scale
# q = obj.rotation_quaternion
# x = obj.location

def scale(a, b):
    x = a[0] * b[0]
    y = a[1] * b[1]
    z = a[2] * b[2]
    return Vector((x, y, z))


def dot(a, b):
    return a @ b

def cross(a, b):
    return a.cross(b)

def normalize(a):
    return a.normalized()

def abs(a):
    if a < 0.0:
        return -a
    return a

def signft(a):
    if a < 0.0:
        return -1.0;
    return 1.0

def conjugate(a):
    return a.conjugated()

ex = Vector((1.0, 0.0, 0.0))
ey = Vector((0.0, 1.0, 0.0))
ez = Vector((0.0, 0.0, 1.0))

def semiaxes_local(obj):
    h = [
        scale(obj.scale, Vector((1.0, 0.0, 0.0))),
        scale(obj.scale, Vector((0.0, 1.0, 0.0))),
        scale(obj.scale, Vector((0.0, 0.0, 1.0))),
    ]
    return h

def semiaxes_oriented(obj):
    hl = semiaxes_local(obj)
    h = [(obj.rotation_quaternion @ hl[i]) for i in range(3)]
    return h

@dataclass
class RigidCuboid:
    x: Vector
    s: Vector
    q: Quaternion
    h: List[Vector]
    
@dataclass
class IntRes:
    inter: bool
    plane: bool
    p1: Vector
    p2: Vector
    n: Vector
    e1: Vector
    e2: Vector
    x1: Vector
    x2: Vector


def boxBoxIntersection(c1, c2):    
    a = [ex, ey, ez]

    a1 = c1.h
    a2 = c2.h

    c = c1.x - c2.x;

    minDist = -float("inf")
    #v3 minN;
    #v3 minA;
    #plane = false;
    #c1PenC2 = true;
    e1 = Vector()
    e2 = Vector()

    # Flächen c1
    for i in range(3):
        n = normalize(a1[i])
        dist = abs(dot(n, c))
        
        for k in range(3):
            dist -= abs(dot(n, a1[k]))
            dist -= abs(dot(n, a2[k]))

        if dist > 0.0:
            return IntRes(False, True, Vector(), Vector(), n, Vector(), Vector(), Vector(), Vector())
        elif dist > minDist:
            minDist = dist
            minN = n
            plane = True
            c1_is_cN = True
            cN = c1
            cP = c2
            e1 = c1.x

    # Flächen c2
    for i in range(3):
        n = normalize(a2[i])
        dist = abs(dot(n, c))
        
        for k in range(3):
            dist -= abs(dot(n, a1[k]))
            dist -= abs(dot(n, a2[k]))

        if dist > 0.0:
            return IntRes(False, True, Vector(), Vector(), n, Vector(), Vector(), Vector(), Vector())
        elif dist > minDist:
            minDist = dist
            minN = n
            plane = True
            c1_is_cN = False
            cN = c2
            cP = c1
            e1 = c2.x

    # Kanten
    for i in range(3):
        for j in range(3):
            n = cross(a1[i], a2[j])
            ndotn = dot(n, n)
            if ndotn == 0.0:
                continue
            n = normalize(n)
            dist = abs(dot(n, c))
            
            for k in range(3):
                dist -= abs(dot(n, a1[k]))
                dist -= abs(dot(n, a2[k]))

            if dist > 0.0:
                return IntRes(False, False, Vector(), Vector(), n, a1[i], a2[j], Vector(), Vector())
            elif dist > minDist:
                minDist = dist
                minN = n
                plane = False
                e1 = a1[i]
                e2 = a2[j]

    q1 = Vector((1, 1, 1))
    q2 = Vector((1, 1, 1))
    x1 = Vector()
    x2 = Vector()
    minDist = -minDist
    if plane:
        n = signft(dot(cP.x - cN.x, minN)) * minN;
        
        def sign0(a):
            if a < 0.0:
                return -1.0
            if a > 0.0:
                return 1.0
            return 0.0
        
        dotx = -sign0(dot(n, cP.q @ a[0]))
        doty = -sign0(dot(n, cP.q @ a[1]))
        dotz = -sign0(dot(n, cP.q @ a[2]))
        
        qP = scale(cP.s, Vector((
            dotx,
            doty,
            dotz
        )))
        
        b = cP.x + (cP.q @ qP)
        
        # NOTE:
        # if dot_ == 0.0 then two edge directions are parallel,
        # compute max lower and min upper height in this direction
        # and correct the height if the calculated point in this direction
        
        # TODO: second cube pi/2 rotation faces parallel
        
        def fixdot0(dotp, _h1, _h2):
            nonlocal b, c1, c2
            
            if dotp == 0.0:
                # _h1 and _h2 are parallel, pick any
                xN = normalize(_h1)
                
                # get upper / lower in xN
                max1 = dot(xN, c1.x + _h1)
                min1 = dot(xN, c1.x - _h1)
                
                if min1 < max1:
                    bot1 = min1
                    top1 = max1
                else:
                    bot1 = max1
                    top1 = min1
            
                max2 = dot(xN, c2.x + _h2)
                min2 = dot(xN, c2.x - _h2)
            
                if min2 < max2:
                    bot2 = min2
                    top2 = max2
                else:
                    bot2 = max2
                    top2 = min2
                
                # use avg height
                maxbot = max(bot1, bot2)
                mintop = min(top1, top2)
                mid = 0.5 * (maxbot + mintop)
                
                # sub old height and add new
                b = b - dot(xN, b)*xN + mid*xN
        
        fixdot0(dotx, a1[0], a2[0])
        fixdot0(doty, a1[1], a2[1])
        fixdot0(dotz, a1[2], a2[2])
        
        qP = conjugate(cP.q) @ (b - cP.x);
        qN = conjugate(cN.q) @ (b + minDist*n - cN.x)
        
        if c1_is_cN:
            q1 = qN
            q2 = qP
        else:
            q1 = qP
            q2 = qN
        
    else:
        # TODO: sign0 fix
        
        n = signft(dot(c2.x - c1.x, minN)) * minN;

        x1 = scale(Vector((signft(dot(n, a1[0])), signft(dot(n, a1[1])), signft(dot(n, a1[2])))), c1.s)
        b1 = c1.x + (c1.q @ x1)
        x2 = scale(Vector((-signft(dot(n, a2[0])), -signft(dot(n, a2[1])), -signft(dot(n, a2[2])))), c2.s)
        b2 = c2.x + (c2.q @ x2)
        
        nE = cross(e1, n)
        dotp = dot(nE, e2)
        if dotp == 0.0:
            dotp = 1.0
        p2 = b2 + (dot(nE, b1 - b2)/dotp)*e2
        p1 = p2 + minDist*n

        q1 = conjugate(c1.q) @ (p1 - c1.x);
        q2 = conjugate(c2.q) @ (p2 - c2.x)
    return IntRes(True, plane, q1, q2, n, e1, e2, x1, x2)


@dataclass
class Args:
    t: float
    inter: bool

argsin = Args(time.time(), False)

def draw(args):
    
    # TODO: ...
    ct = time.time()
    dt = ct - args.t
    if dt > 0.1:
        args.t = ct
        recompute = True
    else:
        recompute = False
    
    recompute = True
    
    objects = bpy.data.objects
    o1 = objects[0]
    o2 = objects[1]
    o1.rotation_mode = "QUATERNION"
    o2.rotation_mode = "QUATERNION"
    h1 = semiaxes_oriented(o1)
    h2 = semiaxes_oriented(o2)    
    c1 = RigidCuboid(o1.location, o1.scale, o1.rotation_quaternion, h1)
    c2 = RigidCuboid(o2.location, o2.scale, o2.rotation_quaternion, h2)
    r = boxBoxIntersection(c1, c2)
        
    shader = gpu.shader.from_builtin('3D_SMOOTH_COLOR')
    shader.bind()
    
    vp = []
    vc = []
    
    def addline(p1, p2, c):
        nonlocal vp, vc
        vp += [p1[:], p2[:]]
        vc += [c, c]
    
    if r.inter:
        p1r = c1.q @ r.p1
        p2r = c2.q @ r.p2
        
        p1w = c1.x + p1r
        p2w = c2.x + p2r
        
        addline(c1.x, c2.x, (1.0, 0.0, 0.0, 1.0))
        addline(c1.x, p1w, (0.0, 1.0, 0.0, 1.0))
        addline(c2.x, p2w, (0.0, 0.0, 1.0, 1.0))
        addline(p1w, p2w, (1.0, 1.0, 0.0, 1.0))
        
        if r.plane:
            addline(r.e1, r.e1 + r.n, (1.0, 1.0, 1.0, 1.0))
        else:
            addline(c1.x, c1.x + r.e1, (0.0, 0.5, 0.0, 1.0))
            addline(c2.x, c2.x + r.e2, (0.0, 0.5, 0.0, 1.0))
            addline(c1.x, c1.x + r.n, (0.0, 0.0, 0.0, 1.0))
            addline(c2.x, c2.x + r.n, (0.0, 0.0, 0.0, 1.0))
            addline(c1.x, c1.x + c1.q @ r.x1, (0.0, 0.0, 0.5, 1.0))
            addline(c2.x, c2.x + c2.q @ r.x2, (0.0, 0.0, 0.5, 1.0))
    else:
        addline(c1.x, c2.x, (1.0, 1.0, 1.0, 1.0))
    
    batch = batch_for_shader(shader, 'LINES', {"pos": vp, "color": vc})
    bgl.glLineWidth(3.0)
    batch.draw(shader)

new_handler = bpy.types.SpaceView3D.draw_handler_add(draw, (argsin,), 'WINDOW', 'POST_VIEW')

# https://blender.stackexchange.com/questions/75612/how-do-you-remove-a-draw-handler-after-its-been-added
dns = bpy.app.driver_namespace
if "the_draw_handler" in dns:
    old_handler = dns["the_draw_handler"]
    bpy.types.SpaceView3D.draw_handler_remove(old_handler, 'WINDOW')
    print("handler removed")
else:
    print("no handler removed")

dns["the_draw_handler"] = new_handler