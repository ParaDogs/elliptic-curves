from random import randint

def extended_euclidean_algorithm(a, b):
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = b, a
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s, old_t

def mulinv(n, p):
    gcd, x, y = extended_euclidean_algorithm(n, p)
    assert (n * x + p * y) % p == gcd

    if gcd != 1:
        # Или n равно 0, или p не является простым.
        raise ValueError(
            '{} has no multiplicative inverse '
            'modulo {}'.format(n, p))
    else:
        return x % p

class Elliptic:
    def __init__(self, field_power, a1, a2, a3, a4, a6):
        self.field_power = field_power
        self.a1 = a1 % field_power
        self.a2 = a2 % field_power
        self.a3 = a3 % field_power
        self.a4 = a4 % field_power
        self.a6 = a6 % field_power
        assert (4*a4*a4*a4 + 27*a6*a6) % self.field_power != 0, "Singular elliptic curve"
    
    def contains(self, point):
        x = point.x
        y = point.y
        if x == y == 0: return True
        p = self.field_power
        return (y*y + self.a1*x*y + self.a3*y) % p == (x*x*x + self.a2*x*x + self.a4*x + self.a6) % p

    def __eq__(self, other):
        return self.a1 == other.a1 and self.a2 == other.a2 and self.a3 == other.a3 and self.a4 == other.a4 and self.a6 == other.a6

    def __str__(self):
        return "y^2 + " + str(self.a1) + "xy + " + str(self.a3) + "y = x^3 + " + str(self.a2) + "x^2 + " + str(self.a4) + "x + " + str(self.a6)
    
    def gen_points(self, generator, n):
        assert isinstance(generator, EPoint), "Wrong generator class. It must be EPoint."
        assert generator.curve == self, "The curves should be the same."
        res = []
        g = generator
        for i in range(n):
            res.append(generator)
            generator += g
        return res


class EPoint:
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x
        self.y = y
        assert self.curve.contains(self), "The elliptic curve does not contain this point."

    def __str__(self):
        return '('+ str(self.x)+ ', '+ str(self.y)+ ')'

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.curve == other.curve and self.x == other.x and self.y == other.y

    def negative(self):
        return EPoint(self.curve, self.x, -self.y % self.curve.field_power)

    def __add__(self, other):
        assert isinstance(other, EPoint), "Wrong operand class. It must be EPoint."
        assert self.curve.contains(other), "The curves should be the same."
        p = self.curve.field_power
        a1 = self.curve.a1
        a2 = self.curve.a2
        a3 = self.curve.a3
        a4 = self.curve.a4
        a6 = self.curve.a6
        x1 = self.x
        y1 = self.y
        x2 = other.x
        y2 = other.y
        null = EPoint(self.curve, 0,0)
        if self == null:
            return other
        elif other == null:
            return self
        elif self == other.negative():
            return EPoint(self.curve,0,0)
        elif self != other:
            la = ((y2 - y1)*mulinv(x2 - x1, p)) % p
            nu = ((y1*x2 - y2*x1)*mulinv(x2 - x1, p)) % p
        else:
            la = ((3*x1*x1 + 2*a2*x1 + a4 - a1*y1)*mulinv(2*y1 + a1*x1 + a3, p)) % p
            nu = ((-x1*x1*x1 + a4*x1 + 2*a6 - a3*y1)*mulinv(2*y1 + a1*x1 + a3, p)) % p
        res_x = (la*la + a1*la - a2 - x1 - x2) % p
        res_y = (-(la + a1)*res_x - nu - a3) % p
        return EPoint(self.curve, res_x, res_y)
    
    def __mul__(self, other):
        assert isinstance(other, int), "Wrong operand class. It must be Epoint"
        res = self
        if other == 0: return EPoint(self.curve, 0, 0)
        for i in range(other-1):
            res += self
        return res

    def __rmul__(self, other):
        return self*other

    def __hash__(self):
        return hash(self.curve.field_power*self.x + self.y)

class Elgamal:
    def __init__(self, curve, generator, secret):
        self.curve = curve
        self.generator = generator
        assert isinstance(self.curve, Elliptic), "Curve must be elliptic"
        assert isinstance(self.generator, EPoint), "Generator must be point on elliptic curve"
        assert self.curve.contains(generator), "Curve does not contain this point."
        d = 1
        p = self.generator
        null = EPoint(self.curve,0,0)
        while (p != null):
            p += self.generator
            d += 1
        self.generator_power = d
        self.open_key = secret*self.generator

    def encrypt(self, mes, key):
        assert isinstance(mes, int), "mes must be integer"
        assert isinstance(key, EPoint), "key must be point on elliptic curve"
        assert self.curve.contains(key), "Curve does not contain key."
        k = 523
        # k = randint(1,self.generator_power-1)
        R = k*self.generator
        P = k*key
        e = (mes*P.x) % self.curve.field_power
        return [R, e]

    def decrypt(self, point, enc, secret):
        Q = secret*point
        mes = (enc*mulinv(Q.x, self.curve.field_power)) % self.curve.field_power
        return mes


el = Elliptic(31991,0,0,0,31988,1000)
G = EPoint(el,0,5585)
n = 32089

# generate points
points = el.gen_points(G, n)
for p in points: print(p)

print("Эллиптическая кривая: " + str(el))
print("Генератор: " + str(G))
print("Сгенерировано точек: " + str(n))
print("Уникальных точек (порядок генератора в группе): " + str(len(set(points))))

# Elgamal scheme
c = 5103
eg = Elgamal(el, G, c)
m = 10000
D = EPoint(el,12507,2027)
print(eg.decrypt(*eg.encrypt(m, D), c))
