#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <vector>

using Vec3 = std::array<double, 3>;

constexpr double MU = 398600.44e9;  // Earth's gravitational parameter (m^3/s^2)
constexpr double EARTH_EQUATORIAL_RADIUS = 6378137.0;  // m
constexpr double J2 = 0;
constexpr double DT = 0.1;
constexpr double PI = 3.14159265358979323846;
constexpr double SPACECRAFT_MASS = 1000.0;  // kg
constexpr double THRUST_ACCELERATION = 0.0 / SPACECRAFT_MASS;  // m/s^2

struct OrbitalElements {
    double semi_major_axis;           // a (m)
    double eccentricity;              // e
    double inclination;               // i (rad)
    double right_ascension;           // Omega, RAAN (rad)
    double argument_of_periapsis;     // omega (rad)
    double true_anomaly;              // nu (rad)
};

struct State {
    Vec3 position;
    Vec3 velocity;
};

constexpr double radians(double degrees) {
    return degrees * PI / 180.0;
}

// Edit these six classical orbital elements to set the initial orbit.
const OrbitalElements INITIAL_ELEMENTS = {
    10000000.0,       // semi-major axis a (m)
    0.0,              // eccentricity e
    radians(0.0),    // inclination i
    radians(40.0),    // RAAN Omega
    radians(35.0),    // argument of periapsis omega
    radians(0.0),     // true anomaly nu
};

// One Keplerian period, based on the input semi-major axis.
const double SIM_TIME = 2.0 * PI * std::sqrt(
    INITIAL_ELEMENTS.semi_major_axis *
    INITIAL_ELEMENTS.semi_major_axis *
    INITIAL_ELEMENTS.semi_major_axis / MU
);

struct Options {
    std::string csv_path = "orbit_cpp.csv";
    std::string plot_path = "orbit_plot.html";
    std::string height_plot_path = "height_plot.html";
    bool make_plot = true;
    bool open_plot = true;
};

double norm(const Vec3& value) {
    return std::sqrt(
        value[0] * value[0] +
        value[1] * value[1] +
        value[2] * value[2]
    );
}

State orbital_elements_to_state(const OrbitalElements& elements) {
    const double a = elements.semi_major_axis;
    const double e = elements.eccentricity;
    const double i = elements.inclination;
    const double raan = elements.right_ascension;
    const double arg_periapsis = elements.argument_of_periapsis;
    const double nu = elements.true_anomaly;

    if (a <= 0.0 || e < 0.0 || e >= 1.0) {
        throw std::invalid_argument(
            "This propagator requires a > 0 and 0 <= e < 1."
        );
    }

    // Perifocal-frame state:
    // p = a(1-e^2)
    // r_PQW = p/(1+e cos(nu)) [cos(nu), sin(nu), 0]
    // v_PQW = sqrt(mu/p) [-sin(nu), e+cos(nu), 0]

    const double p = a * (1.0 - e * e);

    const double radius =
        p / (1.0 + e * std::cos(nu));

    const Vec3 r_perifocal = {
        radius * std::cos(nu),
        radius * std::sin(nu),
        0.0
    };

    const double velocity_scale =
        std::sqrt(MU / p);

    const Vec3 v_perifocal = {
        -velocity_scale * std::sin(nu),
        velocity_scale * (e + std::cos(nu)),
        0.0
    };

    // Q = R3(Omega) R1(i) R3(omega)
    // Maps PQW coordinates into ECI coordinates.

    const double cos_raan = std::cos(raan);
    const double sin_raan = std::sin(raan);

    const double cos_i = std::cos(i);
    const double sin_i = std::sin(i);

    const double cos_argp = std::cos(arg_periapsis);
    const double sin_argp = std::sin(arg_periapsis);

    const double q11 =
        cos_raan * cos_argp -
        sin_raan * sin_argp * cos_i;

    const double q12 =
        -cos_raan * sin_argp -
        sin_raan * cos_argp * cos_i;

    const double q21 =
        sin_raan * cos_argp +
        cos_raan * sin_argp * cos_i;

    const double q22 =
        -sin_raan * sin_argp +
        cos_raan * cos_argp * cos_i;

    const double q31 =
        sin_argp * sin_i;

    const double q32 =
        cos_argp * sin_i;

    const auto rotate_to_eci = [=](const Vec3& perifocal) {
        return Vec3{
            q11 * perifocal[0] +
            q12 * perifocal[1],

            q21 * perifocal[0] +
            q22 * perifocal[1],

            q31 * perifocal[0] +
            q32 * perifocal[1]
        };
    };

    return {
        rotate_to_eci(r_perifocal),
        rotate_to_eci(v_perifocal)
    };
}

Vec3 j2_acceleration(const Vec3& position) {
    const double radius = norm(position);

    if (radius == 0.0) {
        return {0.0, 0.0, 0.0};
    }

    const double z_over_r =
        position[2] / radius;

    const double z_over_r_squared =
        z_over_r * z_over_r;

    // J2 acceleration:
    //
    // a_J2 = -(3/2) J2 mu R^2 / r^5 *
    //
    // [
    //   x(1 - 5z^2/r^2)
    //   y(1 - 5z^2/r^2)
    //   z(3 - 5z^2/r^2)
    // ]

    const double factor =
        -1.5 *
        J2 *
        MU *
        EARTH_EQUATORIAL_RADIUS *
        EARTH_EQUATORIAL_RADIUS /
        (
            radius *
            radius *
            radius *
            radius *
            radius
        );

    return {
        factor *
            position[0] *
            (1.0 - 5.0 * z_over_r_squared),

        factor *
            position[1] *
            (1.0 - 5.0 * z_over_r_squared),

        factor *
            position[2] *
            (3.0 - 5.0 * z_over_r_squared)
    };
}

// --------------------------------------------------
// TOTAL ACCELERATION
// --------------------------------------------------

Vec3 total_acceleration(
    const Vec3& position,
    const Vec3& velocity,
    double time
) {
    const double radius = norm(position);

    if (radius == 0.0) {
        throw std::runtime_error(
            "Spacecraft position reached Earth's center."
        );
    }

    // Central gravity:
    //
    // a = -mu r / |r|^3

    const double gravity_scale =
        -MU /
        (
            radius *
            radius *
            radius
        );

    const Vec3 j2 =
        j2_acceleration(position);

    Vec3 thrust = {
        0.0,
        0.0,
        0.0
    };

    // Thrust is opposite the spacecraft velocity
    // from t = 20 seconds through t = 80 seconds.
    //
    // If THRUST_ACCELERATION = 0,
    // thrust remains exactly zero.

    if (
        time >= 20.0 &&
        time <= 80.0 &&
        THRUST_ACCELERATION != 0.0
    ) {
        const double speed =
            norm(velocity);

        if (speed > 0.0) {
            thrust[0] =
                -THRUST_ACCELERATION *
                velocity[0] /
                speed;

            thrust[1] =
                -THRUST_ACCELERATION *
                velocity[1] /
                speed;

            thrust[2] =
                -THRUST_ACCELERATION *
                velocity[2] /
                speed;
        }
    }

    return {
        gravity_scale * position[0] +
        j2[0] +
        thrust[0],

        gravity_scale * position[1] +
        j2[1] +
        thrust[1],

        gravity_scale * position[2] +
        j2[2] +
        thrust[2]
    };
}

// --------------------------------------------------
// RK4 STATE DERIVATIVES
// --------------------------------------------------

struct StateDerivative {
    Vec3 d_position;
    Vec3 d_velocity;
};

StateDerivative state_derivative(
    const State& state,
    double time
) {
    return {
        state.velocity,

        total_acceleration(
            state.position,
            state.velocity,
            time
        )
    };
}

State add_scaled_derivative(
    const State& state,
    const StateDerivative& derivative,
    double scale
) {
    return {
        {
            state.position[0] +
                derivative.d_position[0] * scale,

            state.position[1] +
                derivative.d_position[1] * scale,

            state.position[2] +
                derivative.d_position[2] * scale
        },

        {
            state.velocity[0] +
                derivative.d_velocity[0] * scale,

            state.velocity[1] +
                derivative.d_velocity[1] * scale,

            state.velocity[2] +
                derivative.d_velocity[2] * scale
        }
    };
}

// --------------------------------------------------
// RK4 PROPAGATOR
// --------------------------------------------------

State rk4_step(
    const State& state,
    double time,
    double dt
) {
    // Classical fourth-order Runge-Kutta.
    //
    // State:
    //
    // y = [r, v]
    //
    // dr/dt = v
    // dv/dt = acceleration(r, v, t)
    //
    // Gravity, J2 and thrust are recalculated
    // during every RK4 intermediate stage.

    const StateDerivative k1 =
        state_derivative(
            state,
            time
        );

    const State state2 =
        add_scaled_derivative(
            state,
            k1,
            0.5 * dt
        );

    const StateDerivative k2 =
        state_derivative(
            state2,
            time + 0.5 * dt
        );

    const State state3 =
        add_scaled_derivative(
            state,
            k2,
            0.5 * dt
        );

    const StateDerivative k3 =
        state_derivative(
            state3,
            time + 0.5 * dt
        );

    const State state4 =
        add_scaled_derivative(
            state,
            k3,
            dt
        );

    const StateDerivative k4 =
        state_derivative(
            state4,
            time + dt
        );

    State next = state;

    for (int axis = 0; axis < 3; ++axis) {
        next.position[axis] +=
            dt / 6.0 *
            (
                k1.d_position[axis] +
                2.0 * k2.d_position[axis] +
                2.0 * k3.d_position[axis] +
                k4.d_position[axis]
            );

        next.velocity[axis] +=
            dt / 6.0 *
            (
                k1.d_velocity[axis] +
                2.0 * k2.d_velocity[axis] +
                2.0 * k3.d_velocity[axis] +
                k4.d_velocity[axis]
            );
    }

    return next;
}

// --------------------------------------------------
// SIMULATION
// --------------------------------------------------

std::vector<State> simulate_orbit() {
    const int n_steps =
        static_cast<int>(
            SIM_TIME / DT
        ) + 1;

    std::vector<State> trajectory(
        n_steps
    );

    trajectory[0] =
        orbital_elements_to_state(
            INITIAL_ELEMENTS
        );

    for (
        int i = 0;
        i < n_steps - 1;
        ++i
    ) {
        const double time =
            static_cast<double>(i) *
            DT;

        trajectory[i + 1] =
            rk4_step(
                trajectory[i],
                time,
                DT
            );
    }

    return trajectory;
}

// --------------------------------------------------
// SAVE CSV
// --------------------------------------------------

bool save_trajectory(
    const std::string& path,
    const std::vector<State>& trajectory
) {
    std::ofstream file(path);

    if (!file) {
        return false;
    }

    file
        << "time_s,"
        << "x_m,"
        << "y_m,"
        << "z_m,"
        << "vx_m_s,"
        << "vy_m_s,"
        << "vz_m_s,"
        << "height_m\n";

    file << std::setprecision(10);

    for (
        std::size_t i = 0;
        i < trajectory.size();
        ++i
    ) {
        const double time =
            static_cast<double>(i) *
            DT;

        const Vec3& position =
            trajectory[i].position;

        const Vec3& velocity =
            trajectory[i].velocity;

        const double height =
            norm(position) -
            EARTH_EQUATORIAL_RADIUS;

        file
            << time << ','
            << position[0] << ','
            << position[1] << ','
            << position[2] << ','
            << velocity[0] << ','
            << velocity[1] << ','
            << velocity[2] << ','
            << height << '\n';
    }

    return true;
}

// --------------------------------------------------
// SAVE 3D ORBIT HTML
// --------------------------------------------------

bool save_plot_html(
    const std::string& path,
    const std::vector<State>& trajectory
) {
    std::ofstream file(path);

    if (!file) {
        return false;
    }

    file << R"HTML(
<!doctype html>
<html lang="en">

<head>

<meta charset="utf-8">

<meta
    name="viewport"
    content="width=device-width, initial-scale=1"
>

<title>
Animated 3D Orbit
</title>

<style>

html,
body {
    width: 100%;
    height: 100%;
    margin: 0;
    overflow: hidden;
    background: #ffffff;
    color: #111827;
    font-family: Arial, sans-serif;
}

#plot {
    width: 100vw;
    height: 100vh;
    display: block;
    cursor: grab;
}

#plot:active {
    cursor: grabbing;
}

#controls {
    position: fixed;
    left: 20px;
    top: 20px;
    z-index: 2;

    display: flex;
    gap: 10px;
    align-items: center;

    padding: 10px 12px;

    border: 1px solid #d1d5db;
    border-radius: 12px;

    background:
        rgba(
            255,
            255,
            255,
            0.9
        );

    box-shadow:
        0 8px 24px
        rgba(
            0,
            0,
            0,
            0.12
        );

    user-select: none;
}

#controls button {
    border:
        1px solid
        #9ca3af;

    border-radius: 8px;

    background: #ffffff;

    color: #111827;

    padding:
        6px 10px;

    font:
        14px Arial,
        sans-serif;

    cursor: pointer;
}

#controls button:hover {
    background: #f3f4f6;
}

#controls label,
#controls span {
    color: #374151;

    font:
        13px Arial,
        sans-serif;
}

#speed {
    width: 140px;
}

</style>

</head>

<body>

<canvas id="plot"></canvas>

<div id="controls">

    <button id="playPause">
        Pause
    </button>

    <button id="reset">
        Reset
    </button>

    <label for="speed">
        Speed
    </label>

    <input
        id="speed"
        type="range"
        min="10"
        max="1000"
        value="200"
        step="10"
    >

    <span id="speedLabel">
        200x
    </span>

    <span id="timeLabel">
        t = 0 s
    </span>

</div>

<script>

const trajectory = [
)HTML";

    file << std::setprecision(10);

    for (
        std::size_t i = 0;
        i < trajectory.size();
        ++i
    ) {
        const Vec3& position =
            trajectory[i].position;

        file
            << "    ["
            << position[0]
            << ", "
            << position[1]
            << ", "
            << position[2]
            << "]";

        if (
            i + 1 <
            trajectory.size()
        ) {
            file << ",";
        }

        file << "\n";
    }

    file << R"HTML(
];

const AXIS_LIMIT =
    8000000;

const DT_SECONDS =
    0.1;

const TOTAL_SIM_TIME =
    (trajectory.length - 1) *
    DT_SECONDS;

const canvas =
    document.getElementById(
        "plot"
    );

const ctx =
    canvas.getContext(
        "2d"
    );

const playPauseButton =
    document.getElementById(
        "playPause"
    );

const resetButton =
    document.getElementById(
        "reset"
    );

const speedSlider =
    document.getElementById(
        "speed"
    );

const speedLabel =
    document.getElementById(
        "speedLabel"
    );

const timeLabel =
    document.getElementById(
        "timeLabel"
    );

let yaw = -0.75;
let pitch = 0.55;
let zoom = 1.0;

let dragging = false;
let lastX = 0;
let lastY = 0;

let playing = true;

let animationTime = 0;

let speedMultiplier =
    Number(
        speedSlider.value
    );

let lastFrameTime = null;

// --------------------------------------------------
// RESIZE
// --------------------------------------------------

function resize() {
    const dpr =
        window.devicePixelRatio ||
        1;

    canvas.width =
        Math.floor(
            window.innerWidth *
            dpr
        );

    canvas.height =
        Math.floor(
            window.innerHeight *
            dpr
        );

    canvas.style.width =
        `${window.innerWidth}px`;

    canvas.style.height =
        `${window.innerHeight}px`;

    ctx.setTransform(
        dpr,
        0,
        0,
        dpr,
        0,
        0
    );

    draw();
}

// --------------------------------------------------
// ROTATION
// --------------------------------------------------

function rotate(point) {
    const [
        x,
        y,
        z
    ] = point;

    const cy =
        Math.cos(yaw);

    const sy =
        Math.sin(yaw);

    const cp =
        Math.cos(pitch);

    const sp =
        Math.sin(pitch);

    const x1 =
        x * cy -
        y * sy;

    const y1 =
        x * sy +
        y * cy;

    const y2 =
        y1 * cp -
        z * sp;

    const z2 =
        y1 * sp +
        z * cp;

    return [
        x1,
        y2,
        z2
    ];
}

// --------------------------------------------------
// 3D -> 2D PROJECTION
// --------------------------------------------------

function project(point) {
    const [
        x,
        y,
        z
    ] = rotate(point);

    const distance =
        25000000;

    const perspective =
        distance /
        (
            distance -
            z
        );

    const scale =
        Math.min(
            window.innerWidth,
            window.innerHeight
        ) *
        0.36 *
        zoom /
        AXIS_LIMIT *
        perspective;

    return {
        x:
            window.innerWidth /
            2 +
            x *
            scale,

        y:
            window.innerHeight /
            2 -
            y *
            scale,

        depth:
            z
    };
}

// --------------------------------------------------
// INTERPOLATED POSITION
// --------------------------------------------------

function pointAtTime(
    timeSeconds
) {
    const rawIndex =
        timeSeconds /
        DT_SECONDS;

    const index =
        Math.floor(
            rawIndex
        );

    const nextIndex =
        (index + 1) %
        trajectory.length;

    const amount =
        rawIndex -
        index;

    const a =
        trajectory[index];

    const b =
        trajectory[nextIndex];

    return [
        a[0] +
            (
                b[0] -
                a[0]
            ) *
            amount,

        a[1] +
            (
                b[1] -
                a[1]
            ) *
            amount,

        a[2] +
            (
                b[2] -
                a[2]
            ) *
            amount
    ];
}

// --------------------------------------------------
// DRAW LINE
// --------------------------------------------------

function drawLine(
    a,
    b,
    color,
    width = 1
) {
    const pa =
        project(a);

    const pb =
        project(b);

    ctx.beginPath();

    ctx.moveTo(
        pa.x,
        pa.y
    );

    ctx.lineTo(
        pb.x,
        pb.y
    );

    ctx.strokeStyle =
        color;

    ctx.lineWidth =
        width;

    ctx.stroke();
}

// --------------------------------------------------
// DRAW TEXT
// --------------------------------------------------

function drawText(
    text,
    point,
    color,
    align = "center"
) {
    const projected =
        project(point);

    ctx.fillStyle =
        color;

    ctx.textAlign =
        align;

    ctx.textBaseline =
        "middle";

    ctx.fillText(
        text,
        projected.x,
        projected.y
    );
}

function formatTick(
    value
) {
    return value.toLocaleString(
        "en-US",
        {
            maximumFractionDigits:
                0
        }
    );
}

function formatTime(
    value
) {
    return value.toLocaleString(
        "en-US",
        {
            maximumFractionDigits:
                0
        }
    );
}

// --------------------------------------------------
// AXES
// --------------------------------------------------

function drawAxes() {
    ctx.font =
        "13px Arial";

    drawLine(
        [
            -AXIS_LIMIT,
            0,
            0
        ],
        [
            AXIS_LIMIT,
            0,
            0
        ],
        "#d94848",
        1.5
    );

    drawLine(
        [
            0,
            -AXIS_LIMIT,
            0
        ],
        [
            0,
            AXIS_LIMIT,
            0
        ],
        "#2f9e44",
        1.5
    );

    drawLine(
        [
            0,
            0,
            -AXIS_LIMIT
        ],
        [
            0,
            0,
            AXIS_LIMIT
        ],
        "#1c7ed6",
        1.5
    );

    const ticks = [
        -8000000,
        -4000000,
        0,
        4000000,
        8000000
    ];

    for (
        const tick of ticks
    ) {
        drawLine(
            [
                tick,
                -120000,
                0
            ],
            [
                tick,
                120000,
                0
            ],
            "#d94848",
            1
        );

        drawLine(
            [
                -120000,
                tick,
                0
            ],
            [
                120000,
                tick,
                0
            ],
            "#2f9e44",
            1
        );

        drawLine(
            [
                -120000,
                0,
                tick
            ],
            [
                120000,
                0,
                tick
            ],
            "#1c7ed6",
            1
        );

        drawText(
            formatTick(tick),
            [
                tick,
                -520000,
                0
            ],
            "#8b1e1e"
        );

        drawText(
            formatTick(tick),
            [
                -520000,
                tick,
                0
            ],
            "#1d6b31"
        );

        drawText(
            formatTick(tick),
            [
                -520000,
                0,
                tick
            ],
            "#155fa8"
        );
    }

    ctx.font =
        "16px Arial";

    drawText(
        "X (m)",
        [
            AXIS_LIMIT *
            1.08,
            0,
            0
        ],
        "#8b1e1e"
    );

    drawText(
        "Y (m)",
        [
            0,
            AXIS_LIMIT *
            1.08,
            0
        ],
        "#1d6b31"
    );

    drawText(
        "Z (m)",
        [
            0,
            0,
            AXIS_LIMIT *
            1.08
        ],
        "#155fa8"
    );
}

// --------------------------------------------------
// TRAJECTORY
// --------------------------------------------------

function drawTrajectory() {
    ctx.beginPath();

    for (
        let i = 0;
        i <
        trajectory.length;
        i++
    ) {
        const p =
            project(
                trajectory[i]
            );

        if (i === 0) {
            ctx.moveTo(
                p.x,
                p.y
            );
        } else {
            ctx.lineTo(
                p.x,
                p.y
            );
        }
    }

    ctx.strokeStyle =
        "#2563eb";

    ctx.lineWidth =
        2.5;

    ctx.stroke();

    const currentIndex =
        Math.floor(
            animationTime /
            DT_SECONDS
        );

    ctx.beginPath();

    for (
        let i = 0;
        i <= currentIndex;
        i++
    ) {
        const p =
            project(
                trajectory[i]
            );

        if (i === 0) {
            ctx.moveTo(
                p.x,
                p.y
            );
        } else {
            ctx.lineTo(
                p.x,
                p.y
            );
        }
    }

    ctx.strokeStyle =
        "#f97316";

    ctx.lineWidth =
        4;

    ctx.stroke();

    const start =
        project(
            trajectory[0]
        );

    const end =
        project(
            trajectory[
                trajectory.length -
                1
            ]
        );

    ctx.fillStyle =
        "#16a34a";

    ctx.beginPath();

    ctx.arc(
        start.x,
        start.y,
        5,
        0,
        Math.PI *
        2
    );

    ctx.fill();

    ctx.fillStyle =
        "#dc2626";

    ctx.beginPath();

    ctx.arc(
        end.x,
        end.y,
        5,
        0,
        Math.PI *
        2
    );

    ctx.fill();
}

// --------------------------------------------------
// SATELLITE
// --------------------------------------------------

function drawSatellite() {
    const satellite =
        project(
            pointAtTime(
                animationTime
            )
        );

    ctx.shadowColor =
        "rgba(249, 115, 22, 0.55)";

    ctx.shadowBlur =
        18;

    ctx.fillStyle =
        "#f97316";

    ctx.beginPath();

    ctx.arc(
        satellite.x,
        satellite.y,
        8,
        0,
        Math.PI *
        2
    );

    ctx.fill();

    ctx.shadowBlur =
        0;

    ctx.strokeStyle =
        "#111827";

    ctx.lineWidth =
        2;

    ctx.beginPath();

    ctx.arc(
        satellite.x,
        satellite.y,
        8,
        0,
        Math.PI *
        2
    );

    ctx.stroke();
}

// --------------------------------------------------
// HUD
// --------------------------------------------------

function drawHud() {
    ctx.fillStyle =
        "#111827";

    ctx.font =
        "20px Arial";

    ctx.textAlign =
        "center";

    ctx.textBaseline =
        "top";

    ctx.fillText(
        "3D Orbit of the Satellite",
        window.innerWidth /
        2,
        18
    );

    ctx.font =
        "13px Arial";

    ctx.textAlign =
        "left";

    ctx.fillStyle =
        "#4b5563";

    ctx.fillText(
        "Drag to rotate, scroll to zoom",
        20,
        window.innerHeight -
        30
    );
}

// --------------------------------------------------
// DRAW
// --------------------------------------------------

function draw() {
    ctx.clearRect(
        0,
        0,
        window.innerWidth,
        window.innerHeight
    );

    drawAxes();

    drawTrajectory();

    drawSatellite();

    drawHud();

    timeLabel.textContent =
        `t = ${formatTime(
            animationTime
        )} s`;
}

// --------------------------------------------------
// ANIMATION
// --------------------------------------------------

function animate(
    timestamp
) {
    if (
        lastFrameTime ===
        null
    ) {
        lastFrameTime =
            timestamp;
    }

    const elapsedSeconds =
        (
            timestamp -
            lastFrameTime
        ) /
        1000;

    lastFrameTime =
        timestamp;

    if (playing) {
        animationTime +=
            elapsedSeconds *
            speedMultiplier;

        while (
            animationTime >=
            TOTAL_SIM_TIME
        ) {
            animationTime -=
                TOTAL_SIM_TIME;
        }

        draw();
    }

    requestAnimationFrame(
        animate
    );
}

// --------------------------------------------------
// CONTROLS
// --------------------------------------------------

playPauseButton.addEventListener(
    "click",
    () => {
        playing =
            !playing;

        playPauseButton.textContent =
            playing
                ? "Pause"
                : "Play";
    }
);

resetButton.addEventListener(
    "click",
    () => {
        animationTime =
            0;

        draw();
    }
);

speedSlider.addEventListener(
    "input",
    () => {
        speedMultiplier =
            Number(
                speedSlider.value
            );

        speedLabel.textContent =
            `${speedMultiplier}x`;
    }
);

// --------------------------------------------------
// CAMERA CONTROLS
// --------------------------------------------------

canvas.addEventListener(
    "mousedown",
    event => {
        dragging =
            true;

        lastX =
            event.clientX;

        lastY =
            event.clientY;
    }
);

window.addEventListener(
    "mouseup",
    () => {
        dragging =
            false;
    }
);

window.addEventListener(
    "mousemove",
    event => {
        if (!dragging) {
            return;
        }

        yaw +=
            (
                event.clientX -
                lastX
            ) *
            0.008;

        pitch +=
            (
                event.clientY -
                lastY
            ) *
            0.008;

        pitch =
            Math.max(
                -1.35,
                Math.min(
                    1.35,
                    pitch
                )
            );

        lastX =
            event.clientX;

        lastY =
            event.clientY;

        draw();
    }
);

canvas.addEventListener(
    "wheel",
    event => {
        event.preventDefault();

        zoom *=
            event.deltaY >
            0
                ? 0.9
                : 1.1;

        zoom =
            Math.max(
                0.4,
                Math.min(
                    3.0,
                    zoom
                )
            );

        draw();
    },
    {
        passive: false
    }
);

window.addEventListener(
    "resize",
    resize
);

resize();

requestAnimationFrame(
    animate
);

</script>

</body>

</html>
)HTML";

    return true;
}

// --------------------------------------------------
// SAVE HEIGHT GRAPH
// --------------------------------------------------

bool save_height_plot_html(
    const std::string& path,
    const std::vector<State>& trajectory
) {
    std::ofstream file(path);

    if (!file) {
        return false;
    }

    file << R"HTML(
<!doctype html>

<html lang="en">

<head>

<meta charset="utf-8">

<meta
    name="viewport"
    content="width=device-width, initial-scale=1"
>

<title>
Spacecraft Height vs. Time
</title>

<style>

html,
body {
    margin: 0;
    width: 100%;
    height: 100%;
    font-family:
        Arial,
        sans-serif;
    background: #fff;
}

#chart {
    width: 100vw;
    height: 100vh;
    display: block;
}

</style>

</head>

<body>

<canvas id="chart"></canvas>

<script>

const samples = [
)HTML";

    file <<
        std::setprecision(10);

    for (
        std::size_t i = 0;
        i < trajectory.size();
        ++i
    ) {
        const double time =
            static_cast<double>(i) *
            DT;

        const double height =
            norm(
                trajectory[i].position
            ) -
            EARTH_EQUATORIAL_RADIUS;

        file
            << "["
            << time
            << ","
            << height
            << "]";

        if (
            i + 1 <
            trajectory.size()
        ) {
            file << ",";
        }

        file << '\n';
    }

    file << R"HTML(
];

const canvas =
    document.getElementById(
        "chart"
    );

const ctx =
    canvas.getContext(
        "2d"
    );

function draw() {
    const dpr =
        window.devicePixelRatio ||
        1;

    const width =
        window.innerWidth;

    const height =
        window.innerHeight;

    canvas.width =
        Math.floor(
            width *
            dpr
        );

    canvas.height =
        Math.floor(
            height *
            dpr
        );

    canvas.style.width =
        width +
        "px";

    canvas.style.height =
        height +
        "px";

    ctx.setTransform(
        dpr,
        0,
        0,
        dpr,
        0,
        0
    );

    const margin = {
        left: 90,
        right: 30,
        top: 60,
        bottom: 70
    };

    const plotWidth =
        width -
        margin.left -
        margin.right;

    const plotHeight =
        height -
        margin.top -
        margin.bottom;

    const times =
        samples.map(
            p => p[0]
        );

    const heights =
        samples.map(
            p => p[1]
        );

    const minTime =
        Math.min(
            ...times
        );

    const maxTime =
        Math.max(
            ...times
        );

    const rawMinHeight =
        Math.min(
            ...heights
        );

    const rawMaxHeight =
        Math.max(
            ...heights
        );

    const rawHeightRange =
        rawMaxHeight -
        rawMinHeight;

    let minHeight;
    let maxHeight;

    // --------------------------------------------------
    // IMPORTANT FIX:
    //
    // Do not massively magnify tiny floating-point
    // numerical errors.
    //
    // With:
    //
    // eccentricity = 0
    // J2 = 0
    // thrust = 0
    //
    // the physical altitude should be constant.
    //
    // RK4 makes the actual numerical error extremely
    // small. If the remaining range is below 1 meter,
    // display a fixed 2-meter vertical range instead
    // of stretching microscopic numerical noise over
    // the full screen.
    // --------------------------------------------------

    if (
        rawHeightRange <
        1.0
    ) {
        const centerHeight =
            0.5 *
            (
                rawMinHeight +
                rawMaxHeight
            );

        minHeight =
            centerHeight -
            1.0;

        maxHeight =
            centerHeight +
            1.0;
    } else {
        const padding =
            rawHeightRange *
            0.08;

        minHeight =
            rawMinHeight -
            padding;

        maxHeight =
            rawMaxHeight +
            padding;
    }

    const x =
        t =>
            margin.left +
            (
                t -
                minTime
            ) /
            (
                maxTime -
                minTime
            ) *
            plotWidth;

    const y =
        h =>
            margin.top +
            (
                maxHeight -
                h
            ) /
            (
                maxHeight -
                minHeight
            ) *
            plotHeight;

    ctx.clearRect(
        0,
        0,
        width,
        height
    );

    // --------------------------------------------------
    // TITLE
    // --------------------------------------------------

    ctx.fillStyle =
        "#111827";

    ctx.font =
        "bold 22px Arial";

    ctx.textAlign =
        "center";

    ctx.fillText(
        "Spacecraft Height vs. Time",
        width /
        2,
        32
    );

    // --------------------------------------------------
    // GRID
    // --------------------------------------------------

    ctx.font =
        "12px Arial";

    ctx.strokeStyle =
        "#e5e7eb";

    ctx.fillStyle =
        "#374151";

    for (
        let i = 0;
        i <= 5;
        ++i
    ) {
        const px =
            margin.left +
            i *
            plotWidth /
            5;

        const time =
            minTime +
            i *
            (
                maxTime -
                minTime
            ) /
            5;

        ctx.beginPath();

        ctx.moveTo(
            px,
            margin.top
        );

        ctx.lineTo(
            px,
            margin.top +
            plotHeight
        );

        ctx.stroke();

        ctx.textAlign =
            "center";

        ctx.fillText(
            time.toFixed(0),
            px,
            margin.top +
            plotHeight +
            22
        );

        const py =
            margin.top +
            i *
            plotHeight /
            5;

        const value =
            maxHeight -
            i *
            (
                maxHeight -
                minHeight
            ) /
            5;

        ctx.beginPath();

        ctx.moveTo(
            margin.left,
            py
        );

        ctx.lineTo(
            margin.left +
            plotWidth,
            py
        );

        ctx.stroke();

        ctx.textAlign =
            "right";

        ctx.fillText(
            value.toFixed(0),
            margin.left -
            10,
            py +
            4
        );
    }

    // --------------------------------------------------
    // BORDER
    // --------------------------------------------------

    ctx.strokeStyle =
        "#111827";

    ctx.lineWidth =
        1.5;

    ctx.strokeRect(
        margin.left,
        margin.top,
        plotWidth,
        plotHeight
    );

    // --------------------------------------------------
    // HEIGHT LINE
    // --------------------------------------------------

    ctx.strokeStyle =
        "#2563eb";

    ctx.lineWidth =
        2;

    ctx.beginPath();

    samples.forEach(
        (point, i) => {
            if (i === 0) {
                ctx.moveTo(
                    x(point[0]),
                    y(point[1])
                );
            } else {
                ctx.lineTo(
                    x(point[0]),
                    y(point[1])
                );
            }
        }
    );

    ctx.stroke();

    // --------------------------------------------------
    // X AXIS LABEL
    // --------------------------------------------------

    ctx.fillStyle =
        "#111827";

    ctx.font =
        "15px Arial";

    ctx.textAlign =
        "center";

    ctx.fillText(
        "Time (seconds)",
        margin.left +
        plotWidth /
        2,
        height -
        20
    );

    // --------------------------------------------------
    // Y AXIS LABEL
    // --------------------------------------------------

    ctx.save();

    ctx.translate(
        22,
        margin.top +
        plotHeight /
        2
    );

    ctx.rotate(
        -Math.PI /
        2
    );

    ctx.fillText(
        "Height above Earth (meters)",
        0,
        0
    );

    ctx.restore();
}

window.addEventListener(
    "resize",
    draw
);

draw();

</script>

</body>

</html>
)HTML";

    return true;
}

// --------------------------------------------------
// TEMP OUTPUT PATH
// --------------------------------------------------

std::string temp_output_path() {
#ifdef _MSC_VER

    char* temp =
        nullptr;

    std::size_t temp_size =
        0;

    if (
        _dupenv_s(
            &temp,
            &temp_size,
            "TEMP"
        ) != 0 ||
        temp == nullptr ||
        temp[0] == '\0'
    ) {
        std::free(
            temp
        );

        return "";
    }

    std::string path =
        temp;

    std::free(
        temp
    );

#else

    const char* temp =
        std::getenv(
            "TEMP"
        );

    if (
        temp == nullptr ||
        temp[0] == '\0'
    ) {
        return "";
    }

    std::string path =
        temp;

#endif

    const char last =
        path[
            path.size() -
            1
        ];

    if (
        last != '\\' &&
        last != '/'
    ) {
#ifdef _WIN32
        path += '\\';
#else
        path += '/';
#endif
    }

    return
        path +
        "orbit_cpp.csv";
}

// --------------------------------------------------
// COMMAND LINE OPTIONS
// --------------------------------------------------

Options options_from_args(
    int argc,
    char* argv[]
) {
    Options options;

    for (
        int i = 1;
        i < argc;
        ++i
    ) {
        const std::string arg =
            argv[i];

        if (
            (
                arg == "-o" ||
                arg == "--output"
            ) &&
            i + 1 <
            argc
        ) {
            options.csv_path =
                argv[++i];
        }

        else if (
            arg ==
            "--plot-output" &&
            i + 1 <
            argc
        ) {
            options.plot_path =
                argv[++i];
        }

        else if (
            arg ==
            "--height-plot-output" &&
            i + 1 <
            argc
        ) {
            options.height_plot_path =
                argv[++i];
        }

        else if (
            arg ==
            "--no-plot"
        ) {
            options.make_plot =
                false;

            options.open_plot =
                false;
        }

        else if (
            arg ==
            "--no-open"
        ) {
            options.open_plot =
                false;
        }
    }

    return options;
}

// --------------------------------------------------
// OPEN GENERATED FILE
// --------------------------------------------------

bool open_file(
    const std::string& path
) {
#ifdef _WIN32

    const std::string command =
        "start \"\" \"" +
        path +
        "\"";

#elif __APPLE__

    const std::string command =
        "open \"" +
        path +
        "\"";

#else

    const std::string command =
        "xdg-open \"" +
        path +
        "\"";

#endif

    return
        std::system(
            command.c_str()
        ) == 0;
}

// --------------------------------------------------
// MAIN
// --------------------------------------------------

int main(
    int argc,
    char* argv[]
) {
    const Options options =
        options_from_args(
            argc,
            argv
        );

    const std::vector<State> trajectory =
        simulate_orbit();

    // --------------------------------------------------
    // CSV
    // --------------------------------------------------

    if (
        !save_trajectory(
            options.csv_path,
            trajectory
        )
    ) {
        const std::string fallback_path =
            temp_output_path();

        if (
            !fallback_path.empty() &&
            save_trajectory(
                fallback_path,
                trajectory
            )
        ) {
            std::cout
                << "Could not write trajectory data to "
                << options.csv_path
                << '\n'

                << "Wrote "
                << trajectory.size()
                << " trajectory samples to "
                << fallback_path
                << '\n';
        } else {
            std::cerr
                << "Could not write trajectory data to "
                << options.csv_path
                << " or the temp directory\n";

            return 1;
        }
    } else {
        std::cout
            << "Wrote "
            << trajectory.size()
            << " trajectory samples to "
            << options.csv_path
            << '\n';
    }

    // --------------------------------------------------
    // PLOTS
    // --------------------------------------------------

    if (
        options.make_plot
    ) {
        if (
            !save_plot_html(
                options.plot_path,
                trajectory
            )
        ) {
            std::cerr
                << "Could not write plot to "
                << options.plot_path
                << '\n';

            return 1;
        }

        std::cout
            << "Wrote animated 3D plot to "
            << options.plot_path
            << '\n';

        if (
            !save_height_plot_html(
                options.height_plot_path,
                trajectory
            )
        ) {
            std::cerr
                << "Could not write height plot to "
                << options.height_plot_path
                << '\n';

            return 1;
        }

        std::cout
            << "Wrote height plot to "
            << options.height_plot_path
            << '\n';

        if (
            options.open_plot &&
            !open_file(
                options.plot_path
            )
        ) {
            std::cerr
                << "Could not open "
                << options.plot_path
                << '\n';

            return 1;
        }
    }

    return 0;
}