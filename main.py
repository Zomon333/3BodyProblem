"""The 3 body problem.

Author: Dagan Poulin
Date: 6/30/24
"""
import math
import tkinter
import turtle
import time
from concurrent.futures import ThreadPoolExecutor
from turtle import Turtle, Vec2D, setup, setworldcoordinates, Screen


def cross(a, b):

    return [
        (a[1] * b[2]) - (a[2] * b[1]),
        (a[2] * b[0]) - (a[0] * b[2]),
        (a[0] * b[1]) - (a[1] * b[0])
    ]


class Body:
    """Celestial body for 3 body problem simulation.

    Only to be used in systems of 3 bodies. Cannot perform n-body simulations
    or simulations of less than 3 bodies.
    """

    def __init__(self, body_color, mass, pos, velocity, g):
        """Initialize a new celestial body with the given parameters."""
        self.G = g

        # Set up the body's turtle for drawing
        self.turtle = Turtle()
        self.turtle.speed(999999)  # We want our turtles to move fast.
        self.turtle.up()
        self.turtle.color(body_color)
        self.turtle.pencolor(body_color)
        self.turtle.shape('circle')
        self.turtle.goto(pos)
        self.turtle.down()
        self.turtle.turtlesize(0.5)

        self.dV = 250

        # Set object's initial mass, position, and velocity.
        self.m = mass
        self.pos = pos
        self.v = velocity
        self.a = Vec2D(0, 0)

    def get_h(self):
        # Because we're in 2D but need a 3D vector, we only have angular
        # momentum around the 3D

        scaled_pos = self.pos * (1/1000)
        scaled_v = self.v * (1/1000)

        pos_3 = [
            scaled_pos[0],
            scaled_pos[1],
            0
        ]
        v_3 = [
            scaled_v[0],
            scaled_v[1],
            0
        ]

        angularish_momentum = cross(pos_3, v_3)
        angularish_momentum_len = get_len3(angularish_momentum)

        earth_grav_param = 3.986 * pow(10, 5)

        tmp = []
        for each in cross(v_3, angularish_momentum):
            tmp.append(each / earth_grav_param)

        pos_n = []
        for each in pos_3:
            pos_n.append(each * (1/get_len3(pos_3)))

        ecc_vec = [
            tmp[0] - pos_n[0],
            tmp[1] - pos_n[1],
            tmp[2] - pos_n[2]
        ]

        tmp = ecc_vec
        ecc_len = get_len3(tmp)
        ecc_vec = []
        for each in tmp:
            ecc_vec.append(each * (1/ecc_len))

        peri = pow(angularish_momentum_len, 2) / (earth_grav_param * (1 + ecc_len))
        apo = pow(angularish_momentum_len, 2) / (earth_grav_param * (1 - ecc_len))

        semi_major = sum([peri, apo]) / 2
        semi_minor = semi_major * math.sqrt(1 - pow(ecc_len, 2))

        return [apo, peri, semi_major, semi_minor]

    def get_data(self):
        """Return the mass and position of the body as a list."""
        return [self.m, self.pos]

    def grav(self, body2, body3):
        """Calculate the gravitational acceleration exerted on this body.

        Accepts two parameters; body2 and body3, both of which are lists
        in the format of [mass, position].
        """
        # Calculate gravitational forces exerted by body2 and body3
        term1 = (-1 * self.G * body2[0] *
                 ((self.pos - body2[1]) *
                  (1 / pow(get_len(self.pos - body2[1]), 3))))
        term2 = (-1 * self.G * body3[0] *
                 ((self.pos - body3[1]) *
                  (1 / pow(get_len(self.pos - body3[1]), 3))))
        return term1 + term2

    def draw(self):
        """Draw the current position."""
        self.turtle.goto(self.pos)

    def step(self, body2, body3, d_t):
        """Update the body's position and velocity for a single time step.

        Accepts three parameters; body2 and body3, both of which are lists
        in the format of [mass, position], and d_t which is the time step.
        """
        self.a = self.grav(body2, body3)
        self.v += self.a * d_t
        self.pos += self.v * d_t


def get_len3(pos):
    """Calculate the magnitude of a Vec3."""
    return math.sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2))


def get_len(pos):
    """Calculate the magnitude of a vector."""
    return math.sqrt(pow(pos[0], 2) + pow(pos[1], 2))


def get_dist(pos1, pos2):
    """Calculate the distance between two points."""
    return abs(get_len(pos2 - pos1))


def main():
    """Driver stub function used to run the simulation."""
    # Simulation constants
    gravity = 6.67 * pow(10, -11)
    radius_moon = 3.84 * pow(10, 8)
    resolution = 1000

    # Set up proper zoom on display
    setup(resolution, resolution)
    setworldcoordinates(-radius_moon * 1.35,
                        -radius_moon * 1.35,
                        radius_moon * 1.35,
                        radius_moon * 1.35)

    axes = Turtle()
    axes.color('white')
    axes.speed(pow(10, 10))
    axes.goto(0, 0)
    axes.down()
    axes.goto(-radius_moon * 4, 0)
    axes.goto(radius_moon * 4, 0)
    axes.goto(0, 0)
    axes.goto(0, -radius_moon * 4)
    axes.goto(0, radius_moon * 4)
    axes.up()

    axes.hideturtle()
    axes.goto(-radius_moon * 1.35, radius_moon * 0.9)
    axes.pd()
    axes.write("Controls:\n"
               "Up & Down: Control Throttle\n"
               "W: Burn Prograde\n"
               "S: Burn Retrograde\n"
               "A: Burn Cicular In\n"
               "D: Burn Cicular Out\n"
               "Space: Pause\n"
               "< & >: Increase/Decrease Time Step\n",
               False,
               'left',
               ('Consolas', 11, 'normal'))

    colors = [
        "green",  # Earth
        "grey",  # Moon
        "red"  # Lagrangian
    ]
    starting_masses = [
        5.97 * pow(10, 24),  # Earth
        7.34 * pow(10, 22),  # Moon
        1000 * pow(10, 0),  # Lagrangian
    ]
    starting_positions = [
        Vec2D(0, 0),  # Earth
        Vec2D(1, 0) * 3.84 * pow(10, 8),  # Moon
        Vec2D(1 / 2, math.sqrt(3) / 2) * 3.84 * pow(10, 8)  # Lagrangian
    ]
    starting_velocities = [
        Vec2D(0, 0),  # Earth
        Vec2D(0, 1) * 1022,  # Moon
        Vec2D(-1 * math.sqrt(3) / 2, 1 / 2) * 1022,  # Lagrangian
    ]

    # Create celestial bodies
    bodies = []
    for i in range(3):
        bodies.append(Body(colors[i],
                           starting_masses[i],
                           starting_positions[i],
                           starting_velocities[i],
                           gravity))

    # Adjust turtle sizes to represent celestial bodies
    screen = Screen()
    earth_image = '.\\earth.gif'
    screen.bgcolor('black')
    moon_image = '.\\moon.gif'
    screen.register_shape(earth_image)
    screen.register_shape(moon_image)

    earth = bodies[0]
    # earth.turtle.turtlesize(6378.1 * pow(10, size_scale))
    earth.turtle.shapesize(0.1, 0.1)
    earth.turtle.resizemode('user')
    earth.turtle.shape(earth_image)

    moon = bodies[1]
    # moon.turtle.turtlesize(1740 * pow(10, size_scale))
    moon.turtle.shape(moon_image)

    lagrange = bodies[2]
    lagrange.turtle.shape('triangle')

    # Simulation duration parameters
    steps = 2048  # Number of steps
    sim_time = 0  # Starting time
    stop = 27 * 60 * 60 * 24  # Duration

    class Timing:
        def __init__(self):
            self.d_t = stop / steps  # Time step
            self.multiplier = 1.0
            self.default_d_t = self.d_t
            self.paused = False

        def toggle_pause(self):
            match self.paused:
                case True:
                    self.paused = False
                    self.d_t = self.default_d_t * self.multiplier
                    pass
                case False:
                    self.paused = True
                    self.d_t = 0
                    pass

        def time_up(self):
            """Increase the time step."""
            self.toggle_pause()
            if round(self.multiplier, 2) < 0.1:
                self.multiplier += 0.01
            else:
                self.multiplier += 0.1
            if round(self.multiplier, 2) > 25:
                self.multiplier = 25
            self.d_t = self.default_d_t * self.multiplier
            self.toggle_pause()

        def time_down(self):
            """Decrease the time step."""
            self.toggle_pause()
            if round(self.multiplier, 2) > 0.1:
                self.multiplier -= 0.1
            else:
                self.multiplier -= 0.01
            if round(self.multiplier, 2) < 0.01:
                self.multiplier = 0.01
            self.d_t = self.default_d_t * self.multiplier
            self.toggle_pause()

    delta = Timing()

    def accel(x, y):
        """Accelerate the spacecraft towards the point clicked."""
        lagrange.v += ((lagrange.pos - Vec2D(x, y)) *
                       -1 *
                       delta.d_t *
                       pow(10, -10))

    def up():
        """Accelerate the spacecraft prograde."""
        lagrange.v += lagrange.v * (1 / get_len(lagrange.v)) * lagrange.dV

    def down():
        """Accelerate the spacecraft retrograde."""
        lagrange.v += lagrange.v * (1 / get_len(lagrange.v)) * lagrange.dV * -1

    def left():
        """Accelerate the spacecraft circular in."""
        prograde = lagrange.v * (1 / get_len(lagrange.v)) * lagrange.dV
        x = prograde[0]
        y = prograde[1]

        circular_in = Vec2D(
            ((x * math.cos(1.5708)) - (y * math.sin(1.5708))),
            ((x * math.sin(1.5708)) + (y * math.cos(1.5708)))
        )

        lagrange.v += circular_in

    def right():
        """Accelerate the spacecraft circular out."""
        prograde = lagrange.v * (1 / get_len(lagrange.v)) * lagrange.dV
        x = prograde[0]
        y = prograde[1]

        circular_in = Vec2D(
            ((x * math.cos(1.5708)) - (y * math.sin(1.5708))),
            ((x * math.sin(1.5708)) + (y * math.cos(1.5708)))
        )

        lagrange.v += (circular_in * -1)

    def throttle_up():
        if lagrange.dV < 25:
            if lagrange.dV < 5:
                lagrange.dV += 1
            else:
                lagrange.dV += 5
        else:
            lagrange.dV += 25

    def throttle_down():
        if lagrange.dV <= 25:
            if lagrange.dV <= 5:
                lagrange.dV -= 1
            else:
                lagrange.dV -= 5
        else:
            lagrange.dV -= 25
        if lagrange.dV < 0:
            lagrange.dV = 0


    turtle.Screen().listen()
    turtle.Screen().onclick(accel)
    turtle.Screen().onkey(up, 'w')
    turtle.Screen().onkey(down, 's')
    turtle.Screen().onkey(left, 'a')
    turtle.Screen().onkey(right, 'd')
    turtle.Screen().onkey(throttle_up, 'Up')
    turtle.Screen().onkey(throttle_down, 'Down')

    turtle.Screen().onkey(delta.toggle_pause, 'space')
    turtle.Screen().onkey(delta.time_up, '.')
    turtle.Screen().onkey(delta.time_down, ',')

    apo_turtle = Turtle()
    apo_turtle.speed(99999999999999999)
    apo_turtle.pu()
    apo_turtle.color('white')
    apo_turtle.goto((radius_moon * 1.1), (-radius_moon) * 1.35)
    apo_turtle.hideturtle()

    ui_turtle = Turtle()
    ui_turtle.speed(999999999999999999)
    ui_turtle.penup()
    ui_turtle.color('white')
    ui_turtle.goto((-radius_moon) * 1.35, (-radius_moon) * 1.35)
    ui_turtle.hideturtle()

    turtle.Screen().tracer(0, 0)

    def perform_offsets(datum):
        """Multithread offset the bodies to be earth-centric."""
        with ThreadPoolExecutor(max_workers=3) as offset_executor:
            def offset(to_offset, pos):
                to_offset.pos -= pos

            offset_executor.submit(offset, earth, datum)
            offset_executor.submit(offset, moon, datum)
            offset_executor.submit(offset, lagrange, datum)

    def perform_steps():
        """Multithread 3 body gravity step."""
        with ThreadPoolExecutor(max_workers=3) as executor:
            executor.submit(earth.step,
                            moon_data,
                            lagrange_data,
                            delta.d_t)
            executor.submit(moon.step,
                            earth_data,
                            lagrange_data,
                            delta.d_t)
            executor.submit(lagrange.step,
                            earth_data,
                            moon_data,
                            delta.d_t)

        for body in [earth, moon, lagrange]:
            body.draw()

    last_frame = 0
    earth.turtle.up()

    # Run the simulation
    while True:  # time < stop:
        start = time.time()
        earth_data = earth.get_data()
        moon_data = moon.get_data()
        lagrange_data = lagrange.get_data()

        # Update positions relative to Earth's position
        perform_offsets(earth_data[1])

        # Update the state of each body
        perform_steps()

        sim_time += delta.d_t

        if time.time() - start > 0.1:
            # if not delta.paused:
            #     if delta.multiplier > 0:
            earth.turtle.clear()
            moon.turtle.clear()
            lagrange.turtle.clear()

        ui_turtle.clear()
        apo_turtle.clear()
        if delta.paused:
            ui_turtle.write(f'Time Warp Percent: PAUSED\nCurrent Throttle:'
                            f' {lagrange.dV} m/s\n\n\nEarth Velocity: <'
                            f'{round(earth.v[0], 2)} m/s, '
                            f'{round(earth.v[1], 2)} m/s>\nEarth Speed: '
                            f'{round(get_len(earth.v), 2)} '
                            f'm/s\n\nMoon Velocity'
                            f': <{round(moon.v[0], 2)} m/s, '
                            f'{round(moon.v[1], 2)} m/s>\nMoon Speed: '
                            f'{round(get_len(moon.v), 2)} m/s\nMoon Altitude: '
                            f'{round(get_dist(moon.pos, earth.pos) * (1/1000), 2)}'
                            f' km\n\nLagrange Velocity: <'
                            f'{round(lagrange.v[0], 2)} m/s, '
                            f'{round(lagrange.v[1], 2)} m/s>\nLagrange Speed: '
                            f'{round(get_len(lagrange.v), 2)} m/s\nLagrange '
                            f'Altitude: '
                            f'{round(get_dist(lagrange.pos, earth.pos) * (1/1000), 2)} km',
                            False,
                            'left',
                            ('Consolas', 11, 'normal'))

            try:
                h = lagrange.get_h()
                apo_turtle.write(f'Orbital Parameters:\n'
                                 f'Apoapsis: {h[0] * 1000} m\n'
                                 f'Periapsis: {h[1] * 1000} m\n'
                                 f'Semi-major axis: {h[2] * 1000} m\n'
                                 f'Semi-minor axis: {h[3] * 1000} m\n',
                                 False,
                                 'right',
                                 ('Consolas', 11, 'normal'))
            except:
                apo_turtle.write(f'NO DATA')

        else:
            ui_turtle.write(f'Time Warp Percent: '
                            f'{round(delta.multiplier * 100, 2)}%\n'
                            f'Current Throttle: '
                            f'{lagrange.dV} m/s\n'
                            f'\n\nEarth Velocity: <'
                            f'{round(earth.v[0], 2)} m/s, '
                            f'{round(earth.v[1], 2)} m/s>\nEarth Speed: '
                            f'{round(get_len(earth.v), 2)}'
                            f' m/s\n\nMoon Velocity: <'
                            f'{round(moon.v[0], 2)} m/s, '
                            f'{round(moon.v[1], 2)} m/s>\nMoon Speed: '
                            f'{round(get_len(moon.v), 2)} m/s\nMoon Altitude: '
                            f'{round(get_dist(moon.pos, earth.pos) * (1/1000), 2)}'
                            f' km\n\nLagrange Velocity: <'
                            f'{round(lagrange.v[0], 2)} m/s, '
                            f'{round(lagrange.v[1], 2)} m/s>\nLagrange Speed: '
                            f'{round(get_len(lagrange.v), 2)} '
                            f'm/s\nLagrange Altitude: '
                            f'{round(get_dist(lagrange.pos, earth.pos) * (1/1000), 2)} km',
                            False,
                            'left',
                            ('Consolas', 11, 'normal'))
            try:
                h = lagrange.get_h()
                apo_turtle.write(f'Orbital Parameters:\n'
                                 f'Apoapsis: {round(h[0], 2)} km\n'
                                 f'Periapsis: {round(h[1], 2)} km\n'
                                 f'Semi-major axis: {round(h[2], 2)} km\n'
                                 f'Semi-minor axis: {round(h[3], 2)} km\n',
                                 False,
                                 'right',
                                 ('Consolas', 11, 'normal'))
            except:
                apo_turtle.write(f'NO DATA')

        if sim_time - last_frame > pow(2, 9.9):
            last_frame = sim_time
            turtle.Screen().update()
        elif delta.paused:
            turtle.Screen().update()
        elif round(abs(delta.d_t), 2) < 0.01:
            turtle.Screen().update()
        time.sleep(abs(delta.d_t) * pow(10, -5.1))


if __name__ == "__main__":
    try:
        main()
    except tkinter.TclError:
        print("TCLError")
        pass
