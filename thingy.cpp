#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

using Vec3 = std::array<double, 3>;

constexpr double MU = 398600.44e9;  // Earth's gravitational parameter (m^3/s^2)
constexpr double RADIUS = 6993000.0;
constexpr double CIRCULAR_SPEED = 7546.106713761982;
constexpr double INCLINATION_DEG = 30.0;
constexpr double DT = 10.0;
constexpr double PI = 3.14159265358979323846;
const double SIM_TIME = 2.0 * PI * std::sqrt((RADIUS * RADIUS * RADIUS) / MU);
const double SpaceCraftMass = 1000.0;  // Mass of the spacecraft in kg
const double ThrustForce = 2000.0 / SpaceCraftMass;  // Thrust force in m/s^2 (2000 N divided by mass)

const Vec3 INITIAL_POSITION = {RADIUS, 0.0, 0.0};
const Vec3 INITIAL_VELOCITY = {
    0.0,
    CIRCULAR_SPEED * std::cos(INCLINATION_DEG * PI / 180.0),
    CIRCULAR_SPEED * std::sin(INCLINATION_DEG * PI / 180.0),
};

struct State {
    Vec3 position;
    Vec3 velocity;
};

struct Options {
    std::string csv_path = "orbit_cpp.csv";
    std::string plot_path = "orbit_plot.html";
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

std::vector<State> simulate_orbit() {
    const int n_steps = static_cast<int>(SIM_TIME / DT);
    std::vector<State> trajectory(n_steps);

    trajectory[0] = {INITIAL_POSITION, INITIAL_VELOCITY};

    for (int i = 0; i < n_steps - 1; ++i) {
        const Vec3& position = trajectory[i].position;
        const Vec3& velocity = trajectory[i].velocity;

        const double radius = norm(position);
        const double gravity_scale = -MU / (radius * radius * radius);

        const double time = i* DT;
        Vec3 currentTrust = {0.0, 0.0, 0.0};

        if (time >= 20.0 && time <=25.0){
            currentTrust[0] = -ThrustForce * velocity[0] / norm(velocity);
            currentTrust[1] = -ThrustForce * velocity[1] / norm(velocity);
            currentTrust[2] = -ThrustForce * velocity[2] / norm(velocity);
        }

        Vec3 acceleration = {
            gravity_scale * position[0] + currentTrust[0],
            gravity_scale * position[1] + currentTrust[1],
            gravity_scale * position[2] + currentTrust[2],
        };

        Vec3 next_velocity = {
            velocity[0] + acceleration[0] * DT,
            velocity[1] + acceleration[1] * DT,
            velocity[2] + acceleration[2] * DT,
        };

        Vec3 next_position = {
            position[0] + next_velocity[0] * DT,
            position[1] + next_velocity[1] * DT,
            position[2] + next_velocity[2] * DT,
        };

        trajectory[i + 1] = {next_position, next_velocity};
    }

    return trajectory;
}

bool save_trajectory(const std::string& path, const std::vector<State>& trajectory) {
    std::ofstream file(path);
    if (!file) {
        return false;
    }

    file << "time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s\n";
    file << std::setprecision(10);

    for (std::size_t i = 0; i < trajectory.size(); ++i) {
        const double time = static_cast<double>(i) * DT;
        const Vec3& position = trajectory[i].position;
        const Vec3& velocity = trajectory[i].velocity;

        file
            << time << ','
            << position[0] << ','
            << position[1] << ','
            << position[2] << ','
            << velocity[0] << ','
            << velocity[1] << ','
            << velocity[2] << '\n';
    }

    return true;
}

bool save_plot_html(const std::string& path, const std::vector<State>& trajectory) {
    std::ofstream file(path);
    if (!file) {
        return false;
    }

    file << R"HTML(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>3D Orbit of the Satellite</title>
<style>
html, body {
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
</style>
</head>
<body>
<canvas id="plot"></canvas>
<script>
const trajectory = [
)HTML";

    file << std::setprecision(10);
    for (std::size_t i = 0; i < trajectory.size(); ++i) {
        const Vec3& position = trajectory[i].position;
        file
            << "    ["
            << position[0] << ", "
            << position[1] << ", "
            << position[2] << "]";
        if (i + 1 < trajectory.size()) {
            file << ',';
        }
        file << '\n';
    }

    file << R"HTML(];

const AXIS_LIMIT = 8000000;
const canvas = document.getElementById("plot");
const ctx = canvas.getContext("2d");

let yaw = -0.75;
let pitch = 0.55;
let zoom = 1.0;
let dragging = false;
let lastX = 0;
let lastY = 0;

function resize() {
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.floor(window.innerWidth * dpr);
    canvas.height = Math.floor(window.innerHeight * dpr);
    canvas.style.width = `${window.innerWidth}px`;
    canvas.style.height = `${window.innerHeight}px`;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
}

function rotate(point) {
    const [x, y, z] = point;
    const cy = Math.cos(yaw);
    const sy = Math.sin(yaw);
    const cp = Math.cos(pitch);
    const sp = Math.sin(pitch);

    const x1 = x * cy - y * sy;
    const y1 = x * sy + y * cy;
    const y2 = y1 * cp - z * sp;
    const z2 = y1 * sp + z * cp;

    return [x1, y2, z2];
}

function project(point) {
    const [x, y, z] = rotate(point);

    const distance = 25000000;
    const perspective = distance / (distance - z);

    const scale = Math.min(window.innerWidth, window.innerHeight)
        * 0.36 * zoom / AXIS_LIMIT * perspective;

    return {
        x: window.innerWidth / 2 + x * scale,
        y: window.innerHeight / 2 - y * scale,
    };
}

function drawLine(a, b, color, width = 1) {
    const pa = project(a);
    const pb = project(b);
    ctx.beginPath();
    ctx.moveTo(pa.x, pa.y);
    ctx.lineTo(pb.x, pb.y);
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.stroke();
}

function drawText(text, point, color, align = "center") {
    const projected = project(point);
    ctx.fillStyle = color;
    ctx.textAlign = align;
    ctx.textBaseline = "middle";
    ctx.fillText(text, projected.x, projected.y);
}

function formatTick(value) {
    return value.toLocaleString("en-US", { maximumFractionDigits: 0 });
}

function drawAxes() {
    ctx.font = "13px Arial";

    drawLine([-AXIS_LIMIT, 0, 0], [AXIS_LIMIT, 0, 0], "#d94848", 1.5);
    drawLine([0, -AXIS_LIMIT, 0], [0, AXIS_LIMIT, 0], "#2f9e44", 1.5);
    drawLine([0, 0, -AXIS_LIMIT], [0, 0, AXIS_LIMIT], "#1c7ed6", 1.5);

    const ticks = [-8000000, -4000000, 0, 4000000, 8000000];
    for (const tick of ticks) {
        drawLine([tick, -120000, 0], [tick, 120000, 0], "#d94848", 1);
        drawLine([-120000, tick, 0], [120000, tick, 0], "#2f9e44", 1);
        drawLine([-120000, 0, tick], [120000, 0, tick], "#1c7ed6", 1);

        drawText(formatTick(tick), [tick, -520000, 0], "#8b1e1e");
        drawText(formatTick(tick), [-520000, tick, 0], "#1d6b31");
        drawText(formatTick(tick), [-520000, 0, tick], "#155fa8");
    }

    ctx.font = "16px Arial";
    drawText("X (m)", [AXIS_LIMIT * 1.08, 0, 0], "#8b1e1e");
    drawText("Y (m)", [0, AXIS_LIMIT * 1.08, 0], "#1d6b31");
    drawText("Z (m)", [0, 0, AXIS_LIMIT * 1.08], "#155fa8");
}

function drawTrajectory() {
    ctx.beginPath();
    for (let i = 0; i < trajectory.length; i++) {
        const p = project(trajectory[i]);
        if (i === 0) {
            ctx.moveTo(p.x, p.y);
        } else {
            ctx.lineTo(p.x, p.y);
        }
    }
    ctx.strokeStyle = "#2563eb";
    ctx.lineWidth = 2.5;
    ctx.stroke();

    const start = project(trajectory[0]);
    const end = project(trajectory[trajectory.length - 1]);

    ctx.fillStyle = "#16a34a";
    ctx.beginPath();
    ctx.arc(start.x, start.y, 5, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#dc2626";
    ctx.beginPath();
    ctx.arc(end.x, end.y, 5, 0, Math.PI * 2);
    ctx.fill();
}

function draw() {
    ctx.clearRect(0, 0, window.innerWidth, window.innerHeight);

    ctx.fillStyle = "#111827";
    ctx.font = "20px Arial";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText("3D Orbit of the Satellite", window.innerWidth / 2, 18);

    ctx.font = "13px Arial";
    ctx.textAlign = "left";
    ctx.fillStyle = "#4b5563";
    ctx.fillText("Drag to rotate, scroll to zoom", 20, window.innerHeight - 30);

    drawAxes();
    drawTrajectory();
}

canvas.addEventListener("mousedown", (event) => {
    dragging = true;
    lastX = event.clientX;
    lastY = event.clientY;
});

window.addEventListener("mouseup", () => {
    dragging = false;
});

window.addEventListener("mousemove", (event) => {
    if (!dragging) {
        return;
    }

    yaw += (event.clientX - lastX) * 0.008;
    pitch += (event.clientY - lastY) * 0.008;
    pitch = Math.max(-1.35, Math.min(1.35, pitch));
    lastX = event.clientX;
    lastY = event.clientY;
    draw();
});

canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    zoom *= event.deltaY > 0 ? 0.9 : 1.1;
    zoom = Math.max(0.4, Math.min(3.0, zoom));
    draw();
}, { passive: false });

window.addEventListener("resize", resize);
resize();
</script>
</body>
</html>
)HTML";

    return true;
}

std::string temp_output_path() {
#ifdef _MSC_VER
    char* temp = nullptr;
    std::size_t temp_size = 0;
    if (_dupenv_s(&temp, &temp_size, "TEMP") != 0 || temp == nullptr || temp[0] == '\0') {
        std::free(temp);
        return "";
    }

    std::string path = temp;
    std::free(temp);
#else
    const char* temp = std::getenv("TEMP");
    if (temp == nullptr || temp[0] == '\0') {
        return "";
    }

    std::string path = temp;
#endif
    const char last = path[path.size() - 1];
    if (last != '\\' && last != '/') {
        path += '\\';
    }

    return path + "orbit_cpp.csv";
}

Options options_from_args(int argc, char* argv[]) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            options.csv_path = argv[++i];
        } else if (arg == "--plot-output" && i + 1 < argc) {
            options.plot_path = argv[++i];
        } else if (arg == "--no-plot") {
            options.make_plot = false;
            options.open_plot = false;
        } else if (arg == "--no-open") {
            options.open_plot = false;
        }
    }

    return options;
}

bool open_file(const std::string& path) {
#ifdef _WIN32
    const std::string command = "start \"\" \"" + path + "\"";
#elif __APPLE__
    const std::string command = "open \"" + path + "\"";
#else
    const std::string command = "xdg-open \"" + path + "\"";
#endif
    return std::system(command.c_str()) == 0;
}

int main(int argc, char* argv[]) {
    const Options options = options_from_args(argc, argv);
    const std::vector<State> trajectory = simulate_orbit();

    if (!save_trajectory(options.csv_path, trajectory)) {
        const std::string fallback_path = temp_output_path();
        if (!fallback_path.empty() && save_trajectory(fallback_path, trajectory)) {
            std::cout
                << "Could not write trajectory data to " << options.csv_path << '\n'
                << "Wrote " << trajectory.size()
                << " trajectory samples to " << fallback_path << '\n';
        } else {
            std::cerr
                << "Could not write trajectory data to " << options.csv_path
                << " or the temp directory\n";
            return 1;
        }
    } else {
        std::cout
            << "Wrote " << trajectory.size()
            << " trajectory samples to " << options.csv_path << '\n';
    }

    if (options.make_plot) {
        if (!save_plot_html(options.plot_path, trajectory)) {
            std::cerr << "Could not write plot to " << options.plot_path << '\n';
            return 1;
        }

        std::cout << "Wrote 3D plot to " << options.plot_path << '\n';

        if (options.open_plot && !open_file(options.plot_path)) {
            std::cerr << "Could not open " << options.plot_path << '\n';
            return 1;
        }
    }

    return 0;
}
