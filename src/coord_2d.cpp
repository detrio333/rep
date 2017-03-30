#include "../headers/coord_2d.h"
#include "math.h"

HAC_USING_LIBRARY_NAMESPACE;

float coord_2d::operator * (coord_2d const& coord2)const
{
	return (this->r*coord2.r + this->z*coord2.z);
};

coord_2d coord_2d::operator+(coord_2d const &__vl)const
{
	coord_2d tmp;
	tmp.r = r + __vl.r;
	tmp.z = z + __vl.z;

	return tmp;
};

coord_2d coord_2d::operator-(coord_2d const &__vl)const
{
	coord_2d tmp;
	tmp.r = r - __vl.r;
	tmp.z = z - __vl.z;

	return tmp;
};

coord_2d coord_2d::operator* (float const &__vl)const
{
	coord_2d tmp;
	tmp.r = r * __vl;
	tmp.z = z * __vl;

	return tmp;
};

coord_2d::coord_2d(float _r, float _z):r(_r), z(_z)
{

};

coord_2d::coord_2d()
{
};

coord_2d coord_2d::operator=(coord_2d &&__vl)
{
	r = __vl.r;
	z = __vl.z;

	return *this;
};

coord_2d coord_2d::operator / (float const &__vl)const
{
	coord_2d tmp;
	tmp.r = r / __vl;
	tmp.z = z / __vl;

	return tmp;
};

float coord_2d::abs()const
{
	return sqrt(this->r*this->r + this->z*this->z);
};
