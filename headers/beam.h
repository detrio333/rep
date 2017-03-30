#ifndef HAC_BEAM_H
#define HAC_BEAM_H

#include "main.h"
#include "bands.h"

#include "list"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT  beam	//! характеристики излучамых лучей
	{
		beam();	//!< констуктор с установкой по умолчанию значений от \f$ -\pi/2 \f$ до \f$ \pi/2 \f$

		bool is_auto;	//!< число лучей, лог (если стоит true, то определет автоматически)
		std::list<float> pow_in_angle;	//<! угол в рад

		space_band_v angle;	//!< угловой раствор излучателя	
	};
}

#endif /*HAC_BEAM_H*/
