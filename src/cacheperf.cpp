#include <intrin.h>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>

namespace {

struct Pos3 {
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
};

struct Vel3 {
	float dx = 0.0f;
	float dy = 0.0f;
	float dz = 0.0f;
};

struct Acc3 {
	float ddx = 0.0f;
	float ddy = 0.0f;
	float ddz = 0.0f;
};

struct Pos4 {
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	float w = 0.0f;
};

struct Vel4 {
	float dx = 0.0f;
	float dy = 0.0f;
	float dz = 0.0f;
	float dw = 0.0f;
};

struct Acc4 {
	float ddx = 0.0f;
	float ddy = 0.0f;
	float ddz = 0.0f;
	float ddw = 0.0f;
};

constexpr unsigned const numParticles = 16 * 1024 * 1024;

void
func0()
{
	std::printf("func0\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<Pos3> _pos(numParticles);
	for (auto & pos : _pos) {
		pos.x = (float)(distribution_0_1000000(random_engine));
		pos.y = (float)(distribution_0_1000000(random_engine));
		pos.z = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<Vel3> _vel(numParticles);
	for (auto & vel : _vel) {
		vel.dx = (float)(distribution_n10_10(random_engine));
		vel.dy = (float)(distribution_n10_10(random_engine));
		vel.dz = (float)(distribution_n10_10(random_engine));
	}

	std::vector<Acc3> _acc(numParticles);
	for (auto& acc : _acc) {
		acc.ddx = (float)(distribution_n1_1(random_engine));
		acc.ddy = (float)(distribution_n1_1(random_engine));
		acc.ddz = (float)(distribution_n1_1(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0/60.0);
	for (unsigned i = 0; i < _pos.size(); ++i) {
		_pos[i].x = _pos[i].x + _vel[i].dx * dt;
		_pos[i].y = _pos[i].y + _vel[i].dy * dt;
		_pos[i].z = _pos[i].z + _vel[i].dz * dt;

		_vel[i].dx = _vel[i].dx + _acc[i].ddx * dt;
		_vel[i].dy = _vel[i].dy + _acc[i].ddz * dt;
		_vel[i].dz = _vel[i].dz + _acc[i].ddy * dt;
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func04()
{
	std::printf("func0\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<Pos4> _pos(numParticles);
	for (auto& pos : _pos) {
		pos.x = (float)(distribution_0_1000000(random_engine));
		pos.y = (float)(distribution_0_1000000(random_engine));
		pos.z = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<Vel4> _vel(numParticles);
	for (auto& vel : _vel) {
		vel.dx = (float)(distribution_n10_10(random_engine));
		vel.dy = (float)(distribution_n10_10(random_engine));
		vel.dz = (float)(distribution_n10_10(random_engine));
	}

	std::vector<Acc4> _acc(numParticles);
	for (auto& acc : _acc) {
		acc.ddx = (float)(distribution_n1_1(random_engine));
		acc.ddy = (float)(distribution_n1_1(random_engine));
		acc.ddz = (float)(distribution_n1_1(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < _pos.size(); ++i) {
		_pos[i].x = _pos[i].x + _vel[i].dx * dt;
		_pos[i].y = _pos[i].y + _vel[i].dy * dt;
		_pos[i].z = _pos[i].z + _vel[i].dz * dt;

		_vel[i].dx = _vel[i].dx + _acc[i].ddx * dt;
		_vel[i].dy = _vel[i].dy + _acc[i].ddz * dt;
		_vel[i].dz = _vel[i].dz + _acc[i].ddy * dt;
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}


void
func1()
{
	std::printf("func1\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] += _vx[i] * dt;
		_py[i] += _vy[i] * dt;
		_pz[i] += _vz[i] * dt;

		_vx[i] += _ax[i] * dt;
		_vy[i] += _ay[i] * dt;
		_vz[i] += _az[i] * dt;
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func2()
{
	std::printf("func2\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 4) {
		__m128 const pxV = _mm_load_ps(&_px[i]);
		__m128 const pyV = _mm_load_ps(&_py[i]);
		__m128 const pzV = _mm_load_ps(&_pz[i]);

		__m128 const vxV = _mm_load_ps(&_vx[i]);
		__m128 const vyV = _mm_load_ps(&_vy[i]);
		__m128 const vzV = _mm_load_ps(&_vz[i]);

		__m128 const dtV = _mm_set1_ps(dt);

		__m128 const vxdtV = _mm_mul_ps(vxV, dtV);
		__m128 const vydtV = _mm_mul_ps(vyV, dtV);
		__m128 const vzdtV = _mm_mul_ps(vzV, dtV);

		__m128 const pxV2 = _mm_add_ps(pxV, vxdtV);
		__m128 const pyV2 = _mm_add_ps(pyV, vydtV);
		__m128 const pzV2 = _mm_add_ps(pzV, vzdtV);

		_mm_store_ps(&_px[i], pxV2);
		_mm_store_ps(&_py[i], pyV2);
		_mm_store_ps(&_pz[i], pzV2);

		__m128 const axV = _mm_load_ps(&_ax[i]);
		__m128 const ayV = _mm_load_ps(&_ay[i]);
		__m128 const azV = _mm_load_ps(&_az[i]);

		__m128 const vxV2 = _mm_add_ps(vxV, axV);
		__m128 const vyV2 = _mm_add_ps(vyV, ayV);
		__m128 const vzV2 = _mm_add_ps(vzV, azV);

		_mm_store_ps(&_vx[i], vxV2);
		_mm_store_ps(&_vy[i], vyV2);
		_mm_store_ps(&_vz[i], vzV2);
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func2ptr()
{
	std::printf("func2ptr\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float const* srcpx = _px.data();
	float const* srcpy = _py.data();
	float const* srcpz = _pz.data();

	float * dstpx = _px.data();
	float * dstpy = _py.data();
	float * dstpz = _pz.data();

	float const* srcvx = _vx.data();
	float const* srcvy = _vy.data();
	float const* srcvz = _vz.data();

	float* dstvx = _vx.data();
	float* dstvy = _vy.data();
	float* dstvz = _vz.data();

	float const* srcax = _ax.data();
	float const* srcay = _ay.data();
	float const* srcaz = _az.data();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 4) {
		__m128 const pxV = _mm_load_ps(srcpx);
		__m128 const pyV = _mm_load_ps(srcpy);
		__m128 const pzV = _mm_load_ps(srcpz);

		__m128 const vxV = _mm_load_ps(srcvx);
		__m128 const vyV = _mm_load_ps(srcvy);
		__m128 const vzV = _mm_load_ps(srcvz);

		__m128 const dtV = _mm_set1_ps(dt);

		__m128 const vxdtV = _mm_mul_ps(vxV, dtV);
		__m128 const vydtV = _mm_mul_ps(vyV, dtV);
		__m128 const vzdtV = _mm_mul_ps(vzV, dtV);

		__m128 const pxV2 = _mm_add_ps(pxV, vxdtV);
		__m128 const pyV2 = _mm_add_ps(pyV, vydtV);
		__m128 const pzV2 = _mm_add_ps(pzV, vzdtV);

		_mm_store_ps(dstpx, pxV2);
		_mm_store_ps(dstpy, pyV2);
		_mm_store_ps(dstpz, pzV2);

		__m128 const axV = _mm_load_ps(srcax);
		__m128 const ayV = _mm_load_ps(srcay);
		__m128 const azV = _mm_load_ps(srcaz);

		__m128 const vxV2 = _mm_add_ps(vxV, axV);
		__m128 const vyV2 = _mm_add_ps(vyV, ayV);
		__m128 const vzV2 = _mm_add_ps(vzV, azV);

		_mm_store_ps(dstvx, vxV2);
		_mm_store_ps(dstvy, vyV2);
		_mm_store_ps(dstvz, vzV2);

		srcpx += 4;
		srcpy += 4;
		srcpz += 4;

		dstpx += 4;
		dstpy += 4;
		dstpz += 4;

		srcvx += 4;
		srcvy += 4;
		srcvz += 4;

		dstvx += 4;
		dstvy += 4;
		dstvz += 4;

		srcax += 4;
		srcay += 4;
		srcaz += 4;
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func2tinyloop()
{
	std::printf("func2tinyloop\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 4) {
		__m128 const vxV = _mm_load_ps(&_vx[i]);
		__m128 const vyV = _mm_load_ps(&_vy[i]);
		__m128 const vzV = _mm_load_ps(&_vz[i]);

		__m128 const dtV = _mm_set1_ps(dt);

		__m128 const vxdtV = _mm_mul_ps(vxV, dtV);
		__m128 const vydtV = _mm_mul_ps(vyV, dtV);
		__m128 const vzdtV = _mm_mul_ps(vzV, dtV);

		__m128 const pxV = _mm_load_ps(&_px[i]);
		__m128 const pyV = _mm_load_ps(&_py[i]);
		__m128 const pzV = _mm_load_ps(&_pz[i]);

		__m128 const pxV2 = _mm_add_ps(pxV, vxdtV);
		__m128 const pyV2 = _mm_add_ps(pyV, vydtV);
		__m128 const pzV2 = _mm_add_ps(pzV, vzdtV);

		_mm_store_ps(&_px[i], pxV2);
		_mm_store_ps(&_py[i], pyV2);
		_mm_store_ps(&_pz[i], pzV2);
	}

	for (unsigned i = 0; i < numParticles; i += 4) {
		__m128 const axV = _mm_load_ps(&_ax[i]);
		__m128 const ayV = _mm_load_ps(&_ay[i]);
		__m128 const azV = _mm_load_ps(&_az[i]);

		__m128 const vxV = _mm_load_ps(&_vx[i]);
		__m128 const vyV = _mm_load_ps(&_vy[i]);
		__m128 const vzV = _mm_load_ps(&_vz[i]);

		__m128 const vxV2 = _mm_add_ps(vxV, axV);
		__m128 const vyV2 = _mm_add_ps(vyV, ayV);
		__m128 const vzV2 = _mm_add_ps(vzV, azV);

		_mm_store_ps(&_vx[i], vxV2);
		_mm_store_ps(&_vy[i], vyV2);
		_mm_store_ps(&_vz[i], vzV2);
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func3()
{
	std::printf("func3\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 4) {
		__m128 const dtV = _mm_set1_ps(dt);

		__m128 const vxV = _mm_load_ps(&_vx[i]);
		__m128 const vxdtV = _mm_mul_ps(vxV, dtV);
		__m128 const pxV = _mm_load_ps(&_px[i]);
		__m128 const pxV2 = _mm_add_ps(pxV, vxdtV);
		_mm_store_ps(&_px[i], pxV2);


		__m128 const vyV = _mm_load_ps(&_vy[i]);
		__m128 const vydtV = _mm_mul_ps(vyV, dtV);
		__m128 const pyV = _mm_load_ps(&_py[i]);
		__m128 const pyV2 = _mm_add_ps(pyV, vydtV);
		_mm_store_ps(&_py[i], pyV2);

		__m128 const vzV = _mm_load_ps(&_vz[i]);
		__m128 const vzdtV = _mm_mul_ps(vzV, dtV);
		__m128 const pzV = _mm_load_ps(&_pz[i]);
		__m128 const pzV2 = _mm_add_ps(pzV, vzdtV);
		_mm_store_ps(&_pz[i], pzV2);

		__m128 const axV = _mm_load_ps(&_ax[i]);
		__m128 const vxV2 = _mm_add_ps(vxV, axV);
		_mm_store_ps(&_vx[i], vxV2);

		__m128 const ayV = _mm_load_ps(&_ay[i]);
		__m128 const vyV2 = _mm_add_ps(vyV, ayV);
		_mm_store_ps(&_vy[i], vyV2);


		__m128 const azV = _mm_load_ps(&_az[i]);
		__m128 const vzV2 = _mm_add_ps(vzV, azV);
		_mm_store_ps(&_vz[i], vzV2);
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

void
func4()
{
	std::printf("func4\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 8) {
		__m256 const pxV = _mm256_load_ps(&_px[i]);
		__m256 const pyV = _mm256_load_ps(&_py[i]);
		__m256 const pzV = _mm256_load_ps(&_pz[i]);

		__m256 const vxV = _mm256_load_ps(&_vx[i]);
		__m256 const vyV = _mm256_load_ps(&_vy[i]);
		__m256 const vzV = _mm256_load_ps(&_vz[i]);

		__m256 const dtV = _mm256_set1_ps(dt);

		__m256 const vxdtV = _mm256_mul_ps(vxV, dtV);
		__m256 const vydtV = _mm256_mul_ps(vyV, dtV);
		__m256 const vzdtV = _mm256_mul_ps(vzV, dtV);

		__m256 const pxV2 = _mm256_add_ps(pxV, vxdtV);
		__m256 const pyV2 = _mm256_add_ps(pyV, vydtV);
		__m256 const pzV2 = _mm256_add_ps(pzV, vzdtV);

		_mm256_store_ps(&_px[i], pxV2);
		_mm256_store_ps(&_py[i], pyV2);
		_mm256_store_ps(&_pz[i], pzV2);

		__m256 const axV = _mm256_load_ps(&_ax[i]);
		__m256 const ayV = _mm256_load_ps(&_ay[i]);
		__m256 const azV = _mm256_load_ps(&_az[i]);

		__m256 const vxV2 = _mm256_add_ps(vxV, axV);
		__m256 const vyV2 = _mm256_add_ps(vyV, ayV);
		__m256 const vzV2 = _mm256_add_ps(vzV, azV);

		_mm256_store_ps(&_vx[i], vxV2);
		_mm256_store_ps(&_vy[i], vyV2);
		_mm256_store_ps(&_vz[i], vzV2);
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}


void
func5()
{
	std::printf("func5\n");

	std::random_device random_device;
	std::mt19937 random_engine(random_device());
	std::uniform_int_distribution<long> distribution_0_1000000(0, 1000000);
	std::uniform_int_distribution<long> distribution_n10_10(-10, 10);
	std::uniform_int_distribution<long> distribution_n1_1(-1, 1);

	std::vector<float> _px(numParticles);
	std::vector<float> _py(numParticles);
	std::vector<float> _pz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_px[i] = (float)(distribution_0_1000000(random_engine));
		_py[i] = (float)(distribution_0_1000000(random_engine));
		_pz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _vx(numParticles);
	std::vector<float> _vy(numParticles);
	std::vector<float> _vz(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_vx[i] = (float)(distribution_0_1000000(random_engine));
		_vy[i] = (float)(distribution_0_1000000(random_engine));
		_vz[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::vector<float> _ax(numParticles);
	std::vector<float> _ay(numParticles);
	std::vector<float> _az(numParticles);
	for (unsigned i = 0; i < numParticles; ++i) {
		_ax[i] = (float)(distribution_0_1000000(random_engine));
		_ay[i] = (float)(distribution_0_1000000(random_engine));
		_az[i] = (float)(distribution_0_1000000(random_engine));
	}

	std::chrono::steady_clock::time_point const startTime = std::chrono::steady_clock::now();

	float dt = (float)(1.0 / 60.0);
	for (unsigned i = 0; i < numParticles; i += 8) {
		__m256 const dtV = _mm256_set1_ps(dt);

		__m256 const vxV = _mm256_load_ps(&_vx[i]);
		__m256 const vxdtV = _mm256_mul_ps(vxV, dtV);
		__m256 const pxV = _mm256_load_ps(&_px[i]);
		__m256 const pxV2 = _mm256_add_ps(pxV, vxdtV);
		_mm256_store_ps(&_px[i], pxV2);


		__m256 const vyV = _mm256_load_ps(&_vy[i]);
		__m256 const vydtV = _mm256_mul_ps(vyV, dtV);
		__m256 const pyV = _mm256_load_ps(&_py[i]);
		__m256 const pyV2 = _mm256_add_ps(pyV, vydtV);
		_mm256_store_ps(&_py[i], pyV2);

		__m256 const vzV = _mm256_load_ps(&_vz[i]);
		__m256 const vzdtV = _mm256_mul_ps(vzV, dtV);
		__m256 const pzV = _mm256_load_ps(&_pz[i]);
		__m256 const pzV2 = _mm256_add_ps(pzV, vzdtV);
		_mm256_store_ps(&_pz[i], pzV2);

		__m256 const axV = _mm256_load_ps(&_ax[i]);
		__m256 const vxV2 = _mm256_add_ps(vxV, axV);
		_mm256_store_ps(&_vx[i], vxV2);

		__m256 const ayV = _mm256_load_ps(&_ay[i]);
		__m256 const vyV2 = _mm256_add_ps(vyV, ayV);
		_mm256_store_ps(&_vy[i], vyV2);


		__m256 const azV = _mm256_load_ps(&_az[i]);
		__m256 const vzV2 = _mm256_add_ps(vzV, azV);
		_mm256_store_ps(&_vz[i], vzV2);
	}

	std::chrono::steady_clock::time_point const stopTime = std::chrono::steady_clock::now();
	auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);
	std::printf("Duration: %fms\n", duration.count() / 1000000.0);
}

}

int
main(int argc, char ** argv)
{
	if (argc > 1) {
		func1();
	} else {
		func2tinyloop();
	}

	return 0;
}
