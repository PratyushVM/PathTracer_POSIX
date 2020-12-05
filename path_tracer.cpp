#include <bits/stdc++.h>
#include <stdlib.h>
#include <pthread.h>
using namespace std;

const double INF = 1e20;

// GLOBAL DEFINITION
int h = 768, w = 1024, samps; // height, width and SPP

double erand()
{
    return (double)rand() / RAND_MAX;
}

// structure defining Vector
struct Vec
{
    // components of the vector
    double x, y, z;

    // constructor for Vec
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    // Vector addition
    Vec operator+(const Vec &b) const
    {
        return (Vec(x + b.x, y + b.y, z + b.z));
    }

    // Vector subtraction
    Vec operator-(const Vec &b) const
    {
        return (Vec(x - b.x, y - b.y, z - b.z));
    }

    // Scalar multiplication
    Vec operator*(double b) const
    {
        return (Vec(x * b, y * b, z * b));
    }

    // Vector multiplication
    Vec mult(const Vec &b) const
    {
        return (Vec(x * b.x, y * b.y, z * b.z));
    }

    // Dot Product
    double dot(const Vec &b) const
    {
        return (x * b.x + y * b.y + z * b.z);
    }

    // Normalize Vector
    Vec &norm()
    {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    // Cross product
    Vec operator%(const Vec &b) const
    {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

// structure defining a Ray
struct Ray
{
    Vec o, d; // Ray = o + td -> o,d are constant vectors

    // constructor
    Ray(Vec o_, Vec d_)
    {
        o = o_;
        d = d_;
    }
};

// GLOBAL DECLARATION

// setting up camera position and direction
Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());

Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135;

// Materials used in radiance
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
};

struct Sphere
{
    double rad;  //radius
    Vec p, e, c; //position,emission,color
    Refl_t refl; //reflection type

    // constructor for a sphere
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_)
    {
        rad = rad_;
        p = p_;
        e = e_;
        c = c_;
        refl = refl_;
    }

    // returns distance if intersects, 0 if no intersection
    double intersect(const Ray &r) const
    {
        // t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 solves intersection point of Ray and Sphere
        Vec op = p - r.o;
        double t, eps = 1e-4;                        // eps - epsilon
        double b = op.dot(r.d);                      // 1/2 b from quadratic equation
        double det = b * b - op.dot(op) + rad * rad; // (b^2 - 4ac)/4 : a=1 because ray is normalized

        if (det < 0)
        {
            // ray misses sphere
            return 0;
        }
        else
        {
            // ray hits sphere
            det = sqrt(det);
        }

        if (b - det > eps)
        {
            t = b - det;
        }
        else if (b + det > eps)
        {
            t = b + det;
        }
        else
        {
            t = 0;
        }

        return t;
    }
};

// Scene definition - hardcoded -could be modified to render different images
Sphere spheres[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   //Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), //Right
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         //Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               //Front
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         //Bottom
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),        //Mirror
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),        //Glass
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)     //Light
};
// alternate scene definitions can be found in https://www.kevinbeason.com/smallpt/extraScenes.txt

// clamps value between 0 and 1
inline double clamp(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

// Maps between 0 to 255, takes gamma factor 2.2 into account
inline int toInt(double x)
{
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

// Checks if ray intersects the spheres
inline bool intersect(const Ray &r, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d;
    t = INF;

    for (int i = int(n) - 1; i >= 0; i--)
    {
        d = spheres[i].intersect(r);
        if (d > 0 && d < t)
        {
            t = d;
            id = i;
        }
    }
    return t < INF;
}

Vec radiance(const Ray &r, int depth)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object

    if (!(intersect(r, t, id)))
    {
        return Vec(); // if miss, return black
    }

    const Sphere &obj = spheres[id]; // hit object

    Vec x = r.o + r.d * t, n = (x - obj.p).norm();
    Vec nl;
    if (n.dot(r.d) < 0)
    {
        nl = n;
    }
    else
    {
        nl = n * -1;
    }
    Vec f = obj.c;

    double p = max({f.x, f.y, f.z}); // maximum reflection

    if (++depth > 5)
    {
        if (erand() < p)
        {
            f = f * (1 / p);
        }
        else
        {
            return obj.e; // Russian Roulette
        }
    }

    if (obj.refl == DIFF)
    {
        // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand(), r2 = erand(), r2s = sqrt(r2);
        Vec w = nl, u;

        if (fabs(w.x) > 0.1)
        {
            u = Vec(0, 1) % w;
        }
        else
        {
            u = Vec(1) % w;
        }

        u = u.norm();
        Vec v = w % u;

        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

        return obj.e + f.mult(radiance(Ray(x, d), depth));
    }
    else if (obj.refl == SPEC)
    {
        // Ideal specular reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth));
    }
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric Refraction

    bool into = n.dot(nl) > 0; // Check if ray is going from outside towards inside
    double nc = 1, nt = 1.5, nnt;
    if (into)
    {
        nnt = nc / nt;
    }
    else
    {
        nnt = nt / nc;
    }
    double ddn = r.d.dot(nl), cos2t;

    // Total Internal Reflection TIR
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
    {
        return obj.e + f.mult(radiance(reflRay, depth));
    }

    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();

    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c;
    if (into)
    {
        c = 1 + ddn;
    }
    else
    {
        c = 1 - tdir.dot(n);
    }

    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);

    if (depth > 2)
    {
        if (erand() < P)
        {
            return obj.e + f.mult(radiance(reflRay, depth) * RP);
        }
        else
        {
            return obj.e + f.mult(radiance(Ray(x, tdir), depth) * TP);
        }
    }
    else
    {
        return obj.e + f.mult(radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr);
    }
}

// structure defining argument for appropriate passing to threads
struct args
{
    int id;
    Vec *c;
};

// mutex m;

void *runner(void *arg)
{
    args *item = (args *)arg;
    int y = item->id;
    Vec *cc = item->c;

    Vec r;

    // Loop columns
    for (unsigned short x = 0; x < w; x++)
    {
        // 2x2 subpixel rows
        for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)
        {
            for (int sx = 0; sx < 2; sx++, r = Vec())
            { // 2x2 subpixel cols
                for (int s = 0; s < samps; s++)
                {
                    double r1 = 2 * erand(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                    double r2 = 2 * erand(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                    Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                    r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
                } // Camera rays are pushed forward to start in interior
                cc[i] = cc[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    srand(time(0));

    // samps -> samples
    if (argc == 2)
    {
        samps = atoi(argv[1]) / 4;
    }
    else
    {
        samps = 1;
    }

    Vec *c = new Vec[w * h];

    // setting up arguments for passing to threads
    args A[h];
    for (int i = 0; i < h; i++)
    {
        A[i].id = i;
        A[i].c = c;
    }

    pthread_t threads[h];

    // double time_;
    // time_ = clock();

    // parallel code
    // spawn threads for each row of the image
    for (int i = 0; i < h; i++)
    {
        pthread_create(&threads[i], NULL, runner, &A[i]);
    }

    // join all threads spawned back to master
    for (int i = 0; i < h; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // write image to ppm file
    FILE *f = fopen("Rendered_image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
    {
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
    // time_ = clock() - time_;
    // cout << "Processor Time taken : " << (double)time_ / CLOCKS_PER_SEC << " seconds\n";
}
