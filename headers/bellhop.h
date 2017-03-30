#ifndef HAC_BELLHOP_H
#define HAC_BELLHOP_H

#include "../../GISDataBase/headers/det_hydro.h"

#include "../headers/main.h"
#include "../headers/ray.h"
#include "../headers/beam.h"
#include "../headers/field_gain.h"
#include "../headers/bellhop_hydrology.h"
#include <iostream>

#include <complex>

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT bellhop	//! класс алгоритмов расчета лучевой структуры
	{
		bellhop_hydrology *local_hydrology = NULL;	//!< местная гидрология

		struct calc_param	//! структура для задания условий расчета
		{
			enum class types_coh		//! типы возможных расчетов кооэфициентов затухания
			{
				coherent,		//!< источник когерентного сигнала
				incoherent,		//!< источник некогерентного сигнала
				semicoherent	//!< источник частично когерентного сигнала
			} type_coh = types_coh::coherent;	//!< тип расчета кооэфициентов затухания

			enum class types_attenuation	//! затухание звука
			{
				none,		//!< не считается
				Thorpe,		//!< по формуле Торпа
				database,	//!< значение приходит из базы
			} type_attenuation = types_attenuation::Thorpe;	//!< затухание звука

			enum class types_beam	//! типы лучей
			{
				gauss,	//!< гауссов луч
				geom	//!< геометрический шляповилный 
			} type_beam = types_beam::gauss;	//!< тип луча

			enum class types_source	//! типы источников
			{
				cylindrical,	//!< для цилиндрических координат
				cartesian	//!< для декартовых координат
			} type_source = types_source::cylindrical;	//!< тип луча

			float depth_source;		//!< глубина источника, м
			float dist_max;			//!< максимальная дистанция расчета лучевой структуры, м

			beam src_beam;	//!< описание луче, выходящих из источника
		};

		bellhop();
		~bellhop();

		// расчеты
		error const set_param(calc_param const&);	//!< запрос данных, установка начальных параметров и параметров по умолчанию
		error const calc_ray_field();				//!< функция расчёта лучевой картины поля
		error const calc_clutter_field();			//!< функция расчёта лучевой картины помехи поля
		error const get_level_field(unsigned int &_Nd, unsigned int &_Nr, std::ostream&);	//!< демонстрационная функция расчёта уровня поля для заданной сетки и вывод в поток
		error const get_ray_field(std::ostream&) const;	//!< записывает геомерию лучей в поток
		error const get_ray_field(std::list<std::map<float,float>> &_container) const;	// запись координат лучей в контейнер

		/*!
		\details функция расчёта передаточной характеристики поля для произвольного расположения приёмника (положение источника в расчётах задаётся изначально и неизменно) т.е. выдаётся значение функции Грина для фиксированных параметров
		\param[in] _r горизонтальное расстояние от источника до приёмника, м
		\param[in] _z глубина приёмника (должна быть меньше глубины волновода), м
		\param[out] _field_gain описание поля
		*/
		error const get_gain(
			float& _r	
			, float& _z	
			, field_gain& _field_gain
			);		//!< функция расчёта передаточной характеристики поля для произвольного расположения приёмника

		/*!
		\details демонстрационная функция расчёта передаточной характеристики поля для произвольного расположения приёмника (положение источника в расчётах задаётся изначально и неизменно) т.е. выдаётся значение функции Грина для фиксированных параметров
		\param[in] _r горизонтальное расстояние от источника до приёмника, м
		\param[in] _z глубина приёмника (должна быть меньше глубины волновода), м
		\param[out] _stream поток с описанием лучевой структуры
		*/
		error const get_gain(
			float& _r
			, float& _z
			, std::ostream& _stream
			) const;		//!< функция вывода в поток лучевой структуры

	private:
		/*!
		\details вычисление параметров луча на каждом шаге. Происходит решение диффернициальных уравнений описанных \ref theory
		
		Используется неявный метод Рунге-Кутты второго порядка. 

		Для его реализации на каждом шаге необходимы как минимум две итерации:

		Прогноз:
		\f[
		y'_{n+1}=y_{n}+h f(x_{n},y_{n})
		\f]
		Коррекция:
		\f[
		y'_{n+1}=y_{n}+h \frac{f(x_{n},y_{n})+f(x_{n+1},y'_{n+1})}{2}
		\f]

		\param[in] _angle угол излучения луча из приемника
		*/
		std::vector<step_ray> trace(float _angle);	//!< расчет по отдельному лучу, исходящего в направдение _angle

		/*!
		\details вычисление параметров луча на каждом шаге. Происходит решение диффернициальных уравнений отписанных \ref theory в главе цилиндрические коордианаты
		
		В отличии от trace(float _angle) вычисление траектории просиходит до первого столкновения с поверхностью

		\param[in] _angle угол излучения луча из приемника
		*/
		std::vector<step_ray> clutter_trace(float _angle);	//!< расчет по отдельному лучу, исходящего в направдение _angle

		/*!
		\details вычисление параметров луча на следующем шаге по параметрам предыдущего
		\param[in] _ray параметры луча на предыдущем шаге
		\param[in] _xtop_lower_r заданая точка поверхности, расположенная позади координаты луча
		\param[in] _xtop_upper_r заданая точка поверхности, расположенная за координатой луча
		\param[in] _xbot_lower_r заданая точка дна, расположенная позади координаты луча
		\param[in] _xbot_upper_r заданая точка дна, расположенная за координатой луча
		*/
		step_ray step(const step_ray &_ray,
			const bdryMod::brdy &_xtop_lower_r, const bdryMod::brdy &_xtop_upper_r,
			const bdryMod::brdy &_xbot_lower_r, const bdryMod::brdy &_xbot_upper_r);	//!< отдельный шаг по лучу

		/*!
		\details Уменьшение величины шага сетки по координатам луча для реализации неявного метода Рунге-Кутты и проверки невыхода луча за границы области расчетов
		\param[in] _ray_x координаты луча
		\param[in] _ray_tray направление луча
		\param[in] _c скорость звука в точке
		\param[in] _xtop_lower_r заданая точка поверхности, расположенная позади координаты луча
		\param[in] _xtop_upper_r заданая точка поверхности, расположенная за координатой луча
		\param[in] _xbot_lower_r заданая точка дна, расположенная позади координаты луча
		\param[in] _xbot_upper_r заданая точка дна, расположенная за координатой луча
		\param _h величина шага по сетке
		*/
		void reducestep(const coord_2d &_ray_x, const coord_2d &_ray_tray, const float &_c,
			const bdryMod::brdy &_xtop_lower_r, const bdryMod::brdy &_xtop_upper_r,
			const bdryMod::brdy &_xbot_lower_r, const bdryMod::brdy &_xbot_upper_r,
			float & _h) const;		//!< Уменьшение величины шага сетки по координатам луча 
		
		/*!
		\details Отражение луча от поверхности
		\param[in] _ray параметры луча
		\param[in] _tbdry касательная к поверхности в точке отражения
		\param[in] _nbdry нормаль к поверхности в точке отражения
		\param[in] _kappa кривизна поверхности в точке отражения
		\param[in] _RefCo Набор коэффициентов отражения
		*/
		step_ray reflect_top(const step_ray  _ray, const coord_2d &_tbdry, const coord_2d &_nbdry, const float &_kappa, bdryMod::brdy::RefCoMod &_RefCo);	//!< Отражение луча от поверхности
		/*!
		\details Отражение луча от дна
		\param[in] _ray параметры луча
		\param[in] _tbdry касательная к дну в точке отражения
		\param[in] _nbdry нормаль к дну в точке отражения
		\param[in] _kappa кривизна дна в точке отражения
		\param[in] _RefCo Набор коэффициентов отражения
		*/
		step_ray reflect_bot(const step_ray  _ray, const coord_2d &_tbdry, const coord_2d &_nbdry, const float &_kappa, bdryMod::brdy::RefCoMod &_RefCo);	//!< Отражение луча от дна

		// расчеты уровня по сетке
		void InfluenceGeoHat(const std::vector<step_ray> &_ray, const float &_alpha);
		void InfluenceGeoGaussian(const std::vector<step_ray> &_ray, const float &_alpha);

		// расчеты уровня в точке
		void InfluenceGeoHatOne(std::multimap<float, field_gain::ray_desc> &_rays, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const;
		void InfluenceGeoGaussianOne(std::multimap<float, field_gain::ray_desc> &_rays, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const;

		// собственные лучи в точке
		void InfluenceGeoHatOne(std::ostream& _stream, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const;
		void InfluenceGeoGaussianOne(std::ostream& _stream, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const;

		/*!
		\details Функция пространственного затухания.
		Учитывает потери при поглощении средой и потери на расширение фронта волны

		Теоретическое описание дано в \ref Attenuation

		*/
		void scalep();	//!< Функция пространственного затухания

		coord_2d src_coord;	//!< координата источника
		float rBox, zBox;	//!< границы расчетов
		float deltas;		//!< начальный шаг по глубине
		float omega;		//!< частота умножаенна на 2*Pi
		float Dalpha;
		ss_at_point ss_at_source;		//!< ВРСЗ в точке излучения
		calc_param local_calc_param;	//!< параметры расчетов загнаны в одну структуру
		beam src_beam;	//!< источник
		std::map<float,std::vector<step_ray>> rays_field;	//!< набор лучей распространения звука по дальности на каждый шаг (первый индекс угол излучения)	
		unsigned int NBeamsOpt;	//!< оптимальное число лучей
		float min_step;	//!< минимальный шаг
		float koef_attenuation;	// коэффициент затухания

 		// для демонстрационной функции для множества приемников
		struct SdRdRMod
		{
			int Nrd;	// число глубин приемника
			float *rd;	// глубина приемника

			int Nr;		// Number of ranges
			float *r;
		};
		SdRdRMod SdRdR;
		std::complex<float> **U;

		friend struct clutter;
	};
}

#endif /*HAC_BELLHOP_H*/