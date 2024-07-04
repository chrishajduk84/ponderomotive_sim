class Particle1D:
    def __init__(self, x, v, charge, mass):
        self.x = x
        self.v = v
        self.charge = charge
        self.mass = mass

class Particle2D:
    def __init__(self, x, y, vx, vy, charge, mass):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.charge = charge
        self.mass = mass

class FixedParticle2D:
    def __init__(self, x, y, charge):
        self.x = x
        self.y = y
        self.charge = charge

