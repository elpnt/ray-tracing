#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

struct Vec {
    double x;
    double y;
    double z;

    Vec(double a, double b, double c) : x(a), y(b), z(c) {}

    Vec operator+(Vec other) {
        return Vec(x + other.x, y + other.y, z + other.z);
    }

    Vec operator-(Vec other) {
        return Vec(x - other.x, y - other.y, z - other.z);
    }

    Vec operator*(double a) { return Vec(x * a, y * a, z * a); }

    double dot(Vec other) { return x * other.x + y * other.y + z * other.z; }

    double norm() { return sqrt(dot(*this)); }

    Vec normalize() {
        double a = (*this).norm();
        return Vec(x / a, y / a, z / a);
    }
};

struct Ray {
    Vec o;  // Origin
    Vec d;  // Direction

    Ray(Vec a, Vec b) : o(a), d(b) {}

    Vec get_origin() { return o; }

    Vec get_direction() { return d; }
};

struct Sphere {
    Vec Center;  // Center
    double Radius;
    Vec Color;

    Sphere(Vec cent, double r, Vec col) : Center(cent), Radius(r), Color(col) {}

    Vec get_center() { return Center; }

    Vec get_color() { return Color; }

    Vec get_normal(Vec p) {
        double a = -1 / Radius;
        return Vec((p - Center).x * a, (p - Center).y * a, (p - Center).z * a);
    }

    bool intersect(Ray &ray, double &t) {
        Vec o = ray.get_origin();
        Vec d = ray.get_direction();

        Vec oc = o - Center;

        const double b = 2 * oc.dot(d);
        const double c = oc.dot(oc) - Radius * Radius;
        double delta = b * b - 4 * c;
        if (delta < 1e-4) return false;

        const double t1 = (-b - sqrt(delta)) / 2;
        const double t2 = (-b + sqrt(delta)) / 2;

        t = (t1 < t2) ? t1 : t2;  // get the first intersection only

        return true;
    }
};

struct Tracer {
    vector<Sphere> spheres;

    Tracer() {}

    void add(Sphere sphere) { spheres.push_back(sphere); }

    Vec Shading(Ray ray, Sphere sphere, double &t) {
        Vec color(0, 0, 0);

        if (sphere.intersect(ray, t)) {
            Vec d = ray.get_direction();
            Vec p = ray.get_origin() + d * t;
            Vec n = sphere.get_normal(p);

            p = p.normalize();
            n = n.normalize();

            double facing_ratio = n.dot(d);

            color = sphere.get_color() * (facing_ratio * 0.5);
        }
        return color;
    }

    Vec trace(int x, int y) {
        Vec pix_color(0, 0, 0);

        double min_t = 10000;
        double t;
        Vec ray_origin = Vec(x, y, 0);
        Vec ray_direction = Vec(0, 0, 1);

        Ray ray(ray_origin, ray_direction);
        for (int i = 0; i < spheres.size(); i++) {
            Vec color = Shading(ray, spheres[i], t);
            if (min_t > t) {
                pix_color = color;
                min_t = t;
            }
        }
        return pix_color;
    }
};

int main() {
    // Window size
    const int W = 500;
    const int H = 500;

    Vec pix_col(0, 0, 0);

    Tracer mytracer = Tracer();

    Sphere sphere(Vec(W * 0.5, H * 0.5, 20), 100, Vec(255, 100, 100));

    mytracer.add(sphere);

    ofstream ofs("result.ppm");

    ofs << "P3\n" << W << " " << H << " 255\n";
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            pix_col = mytracer.trace(y, x);
            ofs << (int)pix_col.x << " " << (int)pix_col.y << " "
                << (int)pix_col.z << "\n";
        }
    }

    return 0;
}