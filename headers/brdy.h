#ifndef HAC_BRDY_H
#define HAC_BRDY_H

#include "main.h"
#include "coord_2d.h"

#include "../../GISDataBase/headers/det_hydro.h"

#include <map>

HAC_LIBRARY_NAMESPACE
{
	struct bdryMod	//! структура дна и поверхности
	{
		error const form_curve_ati_brdy(std::map<float, GDB::det::boundary> const &surf_slices);	//!< формирование криволинейной поверхности по заданным точкам
		error const form_curve_ati_brdy(const float &_depth, const float &_rBox, const float &_freq, const int &_heavy_sea);					//!< формирование плоской поверхности по умолчанию
		error const form_curve_bty_brdy(std::map<float, GDB::det::boundary> const &bottom_slices);	//!< формирование криволинейного дна по заданным точкам
		error const form_curve_bty_brdy(const float &_depth, const float &_rBox);					//!< формирование плоского дна по умолчанию

		struct brdy	//! точка в структуре дна и поверхности
		{
			coord_2d x;			//!< глубина или высота в данной точке, м
			coord_2d normal;	//!< нормаль в данной точке
			coord_2d t;			//!< касательная (tangent for a segment)
			float Len;			//!< длина касательной (length of tangent (temporary variable to normalize tangent))
			float Kappa;		//!< кривизна поверхности в данной точке

			coord_2d Nodet;
			coord_2d Noden;

			struct RefCoMod		//! структура коэффициентов отражения
			{
				struct ReflectionCoef	//! коэффициент отражения
				{
					float theta;//!< угло падения, рад
					float R;	//!< отн. ед. (не дБ)
					float phi;	//!< дополнительный фазовый сдвиг
				};

				std::map<float, ReflectionCoef> RefSet;

				ReflectionCoef REFCO(float &_alpha);	//!< возвращает коэффициент отражения по углу
			};
			
			RefCoMod RefCo;	//!< коэффициенты отражения

			std::map<float, float> ScaSet;	//!< коэффициенты рассеяния (угол, рад - коэф, отн. ед)
				//	угол,	коэф, отн.ед
		};

		std::map<float, brdy> map_brdy;	// 1 - расстояние от источника, 2 - ...

		enum class types	//! типы интерполяции
		{
			curve,	//!< криволинейная интерполяция
			linear	//!< кусочнозаданная линейная интерполяция
		};
		types type;	//!< интерполяция
	};
}

#endif /*HAC_BRDY_H*/
