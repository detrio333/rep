#ifndef HAC_COORD_H
#define HAC_COORD_H

#include "main.h"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT coord_2d //! используюется большое количество 2d координат, целесообразно выделить для них свою структуру
	{
		coord_2d();
		coord_2d(float, float);	//!< конструктор с заданием значений

		float r;	//!< горизонтальная координата, м
		float z;	//!< вертикальная координата, м

		float inline operator * (coord_2d const&)const;	//!< скалярное умножение
		inline coord_2d operator = (coord_2d &&);	//!< приравнивание
		inline coord_2d operator + (coord_2d const&) const;	//!< сложение 
		inline coord_2d operator - (coord_2d const&) const;	//!< вычитание
		inline coord_2d operator * (float const&) const;	//!< умножение
		inline coord_2d operator / (float const&) const;	//!< деление
		float inline abs()const;	//!< модуль
	};
};

#endif /*HAC_COORD_H*/