#ifndef HAC_BANDS_H
#define HAC_BANDS_H

#include "main.h"

#include <functional>

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT freq_band	//! частотный диапазон
	{
		freq_band();	//!< конструктор с установкой неопределенных значений

		/*!
		\param __low нижняя граница диапазона, Гц
		\param __high верхняя граница диапазона, Гц
		*/
		freq_band(float const __low, float const __high);	//!< конструктор с установкой значений

		float low;	//!< нижняя граница диапазона, Гц
		float high;	//!< верхняя граница диапазона, Гц

		float const get_freq_avg()const;	//!< возвращает среднегеометрическую частоту диапазона 
	};

	struct HAC_EXPORT time_band	//! временной интервал
	{
		time_band();	//!< конструктор с установкой неопределенных значений

		/*!
		\param __begin начальное время интервала, мс
		\param __end конечное время интервала, мс
		*/
		time_band(unsigned int const __begin, unsigned int const __end);	//!< конструктор с установкой значений

		unsigned int begin;	//!< начальное время интервала, мс
		unsigned int end;	//!< конечное время интервала, мс
	};

	struct HAC_EXPORT space_band_v	//! пространственный диапазон в вертикальной плоскости
	{
		space_band_v();	//!< конструктор с установкой неопределенных значений

		/*!
		\param __low нижняя граница диапазона, рад
		\param __high верхняя граница диапазона, рад
		*/
		space_band_v(float const __low, float const __high);	//!< конструктор с установкой значений

		float low;	//!< нижняя граница диапазона, рад
		float high;	//!< верхняя граница диапазона, рад

		float const get_space_avg()const;	//!< возвращает среднеарифметический угол

		const bool angle_contains(const float &) const;	//!< фунцкия определения принадлежности угла пространственному диаппазону в ВП
		const bool operator == (const space_band_v &) const;	//!< оператор сравнения
	};

	// функтор для хеширования
	struct HAC_EXPORT hash_space_band_v
	{
		size_t operator()(const space_band_v &_space_band_v) const
		{
			return std::hash<float>()(_space_band_v.low) + std::hash<float>()(_space_band_v.high);
		}
	};
};

#endif /*HAC_BANDS_H*/
