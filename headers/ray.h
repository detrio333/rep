#ifndef HAC_RAY_H
#define HAC_RAY_H

#include "main.h"
#include "coord_2d.h"

#include "vector"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT step_ray	//! по аналогии с готовой программой
	{
		float c;		//!< скорость звука в данной точке, м/с
		coord_2d x;		//!< координаты, м
		coord_2d Tray;	//!< направление распространения

		coord_2d p;
		coord_2d q;
		float tau;		//!< delay time (то есть время распространения луча по траектории), с
		float Len;		//!< длина пути, м

		float Rfa;		//!< относительная мощность, ед.
		float phase;	//!< фаза, рад

		enum class events //! решил добавить отметку об отражениях
		{
			surface,
			volume,
			seabed
		};
		events event;
	};
}

#endif /*HAC_RAY_H*/