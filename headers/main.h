#ifndef HAC_MAIN_H
#define HAC_MAIN_H

#define HAC_LIBRARY_NAMESPACE_NAME HAC
#define H HAC_LIBRARY_NAMESPACE_NAME
#define HAC_LIBRARY_NAMESPACE namespace HAC_LIBRARY_NAMESPACE_NAME
#define H HAC_LIBRARY_NAMESPACE_NAME
#define HAC_USING_LIBRARY_NAMESPACE using HAC_LIBRARY_NAMESPACE

#ifdef HAC_DLL
#ifdef HAC_DLL_EXPORT
#define HAC_EXPORT __declspec(dllexport)
#else
#define HAC_EXPORT __declspec(dllimport)
#endif
#else
#define HAC_EXPORT
#endif

HAC_LIBRARY_NAMESPACE //! пространство имен библиотеки расчета лучевой структуры
{
	/*! \namespace HAC_LIBRARY_NAMESPACE_NAME
	пространство имен библиотеки расчета лучевой структуры
	*/
	enum error //! описание ошибок завершения функций
		{
			OK_terminate,				//!< нормальное завершение программы
			field_dont_calculated,		//!< поле не рассчитаное
			brdy_param_not_set,			//!< не заданы параметры дна или поверхности
			spp_param_not_set,			//!< ВРСЗ не задано
			limit_error,				//!< выход за пределы рассчитаного поля
		};

	namespace constants	//! константы расчетов
	{
		unsigned int const MaxSteps = 10000;	//!< максимальное число шагов по х
		unsigned int const MaxNRay = 1000;	//!< максимальное число лучей (пока не используется)
		unsigned int const MaxNFreq = 100;	//!< максимальное число частот в частотном диапазоне (пока не используется)
	}
};

#endif /*HAC_MAIN_H*/


/**
@mainpage Библиотека гидроакустических расчетов

Аннотация

Оглавление:
- \subpage theory
- \subpage data2
- \subpage tr_eq
- \subpage Attenuation

\page theory Теория звукового поля

\section data1 Акустическое волновое уравнение
Рассмотрим волновое уравнение для излучателя в волноводе
\f[
\nabla ^{2} p -  \frac{1}{c^{2}} \frac{\partial ^{2} p}{\partial t^{2}}
\f]

где \f$ p(r,t) \f$ - звуковое давление акустической волны.

После фурье-преобразования получается уравнение Гельмгольца

\f[
(\nabla ^{2} -  \frac{\omega ^{2}}{c^{2}}) P(r,\omega)
\f]

где 

\f[
P(r, \omega) = \int_{-\infty}^{\infty} p(r,t) exp(-i \omega t) d t
\f]

Решение уравнения представляется в виде плоской волны:
\f[
P(\omega) = A exp(-i \omega \tau)
\f]

\f$ А \f$ - медленно изменяющася амплитуда, \f$ \omega \tau \f$ - быстро меняющаяся относительно ампитуды фаза. Поверхность с констатой \f$ \omega \tau \f$ называется волновым фронтом
, поверхность с простоянной \f$ \tau \f$ по аналогии временной фронт.

Подставновка такого решения в уравнение Гельмгольца приводит к лучевой теории в приближении \f$ \frac{\nabla ^{2} A}{A} \ll k^{2} \f$. \f$ k=\frac{\omega}{c} \f$.

Для полученного выражения вещественная и мнимая часть рассматриваются отдельно:

веществаенная (<a href="data2.html">уравнение эйконала</a>, дает геометрию лучей)
\f[
\nabla \tau = \frac{1}{c}
\f]
мнимая (<a href="tr_eq.html">уравнение переноса</a>, дает уменьшение мощности при переносе)
\f[
2 (\nabla A \cdot \nabla \tau) + A \nabla ^{2} \tau = 0
\f]

\page data2 Решение уравнения эйконала

Уравнение эйконала можно переписать в виде

\f[
\frac{d \tau}{ds}=\frac{1}{c}
\f]

где \f$ ds \f$ - дистанция пройденная акустическим лучем.

Для описания распространения волны между точками A и B:
\f[
\tau= \int_{A}^{B} \frac{1}{c} ds
\f]

\section data3 Формализм Лагранжа и принцип Ферма

Перепишем уравнение эйконала согласно формализму Лагранжа

\f[
\tau= \int_{A}^{B} L ds
\f]

Лагранжиан системы является функцией координат и обощенных скоростей:

\f[
L(x,y,z,x',y',z')=\frac{1}{c} \sqrt{x'^{2}+y'^{2}+z'^{2}}
\f]

где \f$ x'=\frac{dx}{ds} \f$, \f$ y'=\frac{dy}{ds} \f$, \f$ z'=\frac{dz}{ds} \f$.

\f{eqnarray*}{
d \tau &=& \ \int_{A}^{B} [(\frac{\partial L}{\partial x}dx+\frac{\partial L}{\partial x'}dx')+(\frac{\partial L}{\partial y}dy+\frac{\partial L}{\partial y'}dy')+(\frac{\partial L}{\partial z}dz+\frac{\partial L}{\partial z'}dz')] ds \\
&=& \ \frac{\partial L}{\partial x'}dx + \frac{\partial L}{\partial y'}dy + \frac{\partial L}{\partial z'}dz  |^B_A + \\
&=& \ \int_{A}^{B} [(\frac{\partial L}{\partial x}dx - \frac{d}{ds}\frac{\partial L}{\partial x'}dx)+(\frac{\partial L}{\partial y}dy-\frac{d}{ds}\frac{\partial L}{\partial y'}dy)+(\frac{\partial L}{\partial z}dz-\frac{d}{ds}\frac{\partial L}{\partial z'}dz)] ds
\f}

Согласно принципу Ферма \f$ d \tau = 0 \f$. Учитывая, что вторая строка уравнения равна 0, получаем

\f[
\frac{d}{ds}\frac{\partial L}{\partial x'} - \frac{\partial L}{\partial x} = 0
\f]
\f[
\frac{d}{ds}\frac{\partial L}{\partial y'} - \frac{\partial L}{\partial y} = 0
\f]
\f[
\frac{d}{ds}\frac{\partial L}{\partial z'} - \frac{\partial L}{\partial z} = 0
\f]

С другой стороны 

\f[
\frac{\partial L}{\partial x'} = \frac{x'}{\sqrt{x'^2+y'^2+z'^2}} L=Lx'=\frac{1}{c} \frac{dx}{ds}
\f]
\f[
\frac{\partial L}{\partial y'} = \frac{y'}{\sqrt{x'^2+y'^2+z'^2}} L=Ly'=\frac{1}{c} \frac{dy}{ds}
\f]
\f[
\frac{\partial L}{\partial z'} = \frac{z'}{\sqrt{x'^2+y'^2+z'^2}} L=Lz'=\frac{1}{c} \frac{dz}{ds}
\f]

При подставновки полученного выражение выше:

\f[
\frac{d}{ds}(\frac{1}{c} \frac{dx}{ds}) = \frac{\partial}{\partial x} (\frac{1}{c})
\f]
\f[
\frac{d}{ds}(\frac{1}{c} \frac{dy}{ds}) = \frac{\partial}{\partial y} (\frac{1}{c})
\f]
\f[
\frac{d}{ds}(\frac{1}{c} \frac{dz}{ds}) = \frac{\partial}{\partial z} (\frac{1}{c})
\f]

Введем величину \f$ \sigma = \frac{1}{c} \f$, с новым параметром 

\f[
\frac{d}{ds}(\sigma \frac{dx}{ds}) = \frac{\partial \sigma}{\partial x}
\f]
\f[
\frac{d}{ds}(\sigma \frac{dy}{ds}) = \frac{\partial \sigma}{\partial y}
\f]
\f[
\frac{d}{ds}(\sigma \frac{dz}{ds}) = \frac{\partial \sigma}{\partial z}
\f]

представление её в векторном виде дает:

\f[
\frac{d \sigma}{ds} = \nabla \sigma
\f]

\section data4 Цилиндрические координаты

При переходе от x,y,z к r,z

\f[
\frac{dr}{ds} = \frac{\sigma_{r}}{\sigma};  \frac{dz}{ds} = \frac{\sigma_{z}}{\sigma}
\f]
\f[
\frac{\sigma_{r}}{ds} = \frac{\partial \sigma}{\partial r};	\frac{\sigma_{z}}{ds} = \frac{\partial \sigma}{\partial z};
\f]

Или в более классическом варианте:

\f[
\frac{dr}{ds} = c(s) \sigma_{r};  \frac{\sigma_{r}}{ds} = - \frac{1}{c^2} \frac{\partial c}{\partial r}
\f]
\f[
\frac{dz}{ds} = c(s) \sigma_{z};  \frac{\sigma_{z}}{ds} = - \frac{1}{c^2} \frac{\partial c}{\partial z}
\f]

Для нахождения траектории движения лучей необходимо найти решения данных уравнений с начальными условиями
\f[
r(0)=r_{0},	z(0)=z_{0},	\sigma_{r}(0)=\frac{cos(\alpha_{0})}{c(0)},	\sigma_{z}(0)=\frac{sin(\alpha_{0})}{c(0)}
\f]

\page tr_eq Решение уравнения переноса

Решение уравнения Эйконала позволяет посчитать Якобиан \f$ J \f$ между декартовыми координатами (x,y,z) и сферическими (s, &alpha;1, &alpha;2)

Рассмотрим единичный вектор \f$ e_{s} \f$ направленнный вдоль луча
\f[
\ e_{s} = \begin{bmatrix} dx/ds \\ dy/ds \\ dz/ds \end{bmatrix}
\f]

Первая часть уравнения переноса может быть переписана как

\f[
\nabla A \cdot \nabla \tau = \frac{dA}{ds}e_{s} \cdot \frac{d \tau}{ds}e_{s} = \frac{dA}{ds}e_{s} \cdot \frac{1}{c}e_{s} = \frac{1}{c} \frac{dA}{ds}
\f]

Для второй части уравнения заметим, что аналитические свойства якобиана позволяют переписать

\f[
\nabla ^{2} \tau = \nabla \cdot \nabla \tau = \frac{1}{J} \frac{d}{ds} \frac{J}{c}
\f]

Общее выражение для уравнения переноса будет выглядеть 

\f[
\frac{2}{c} \frac{dA}{ds} + \frac{A}{J} \frac{d}{ds} \frac{J}{c} = 0
\f]

решение которого 

\f[
A = \frac{A_{0}}{\sqrt{J/c(s)}}
\f]

где \f$ A_{0} \f$ - константа, зависящая от акустического источника.

Решение волнового уравнения может быть записано как

\f[
P(r,\omega) = A_{0} \sqrt{\frac{c(s)}{J}} e^{-i \omega \tau}
\f]

Для сферических координат, применимых к излучению точечного источника
\f[
c(s) \sim c(0) ; J \sim s^{2}cos(\theta (0))
\f]

\f$ \theta (0) \f$ - угол скольжения на выходе из излучателя

С другой стороны, близкое к истонику акустическое поле может быть аппроксимировано сферической волной, для которой справедливо

\f[
A_{0} \sqrt{\frac{c(0)}{s^{2}cos \theta(0)}} e^{-i \omega s/c(0)} = \frac{1}{4 \pi s} e^{-i \omega s/c(0)}
\f]

это ведет к следующему отношению

\f[
A_{0} = \frac{1}{4 \pi} \sqrt{\frac{cos \theta(0)}{c(0)}} 
\f]

Классическое решение волнового уравнения тогда:

\f[
P(r,\omega) = \frac{1}{4 \pi} \sqrt{\frac{c(s)}{с(0)}\frac{cos \theta(0)}{J}} e^{-i \omega \tau}
\f]

К несчастью, классическое решение с использованием якобиана имеет недостатки. Каждый раз проходя точки каустики якобиан обращается в ноль.

\page Attenuation Затухание

\section Att Затухание

Как показано в предыдущих частх теоретического описания общее решение для поля акустического давления:

\f[
P(r,\omega) = A e^{-i \omega \tau}
\f]

Предыдущее выражение не учитывает потери энергии в результате поглощения средой и в результате отражений от границ среды.
Скорректированное значение амплитуды &alpha; будет отличаться от начального значения A на значение пространственного затухания \f$ \phi_{V} \f$
и затухания при отражении \f$ \phi_{r} \f$

\f[
\alpha = A \phi_{V} \phi_{r}
\f]

\subsection boundary Отражение от поверхности

Затухание при отражении от поверхности

Коэффициенты отражения от поверхности задаются из базы на сетке значений

\subsection Volume Пространственное затухание

Пространственное затухание в океане имеет химическую природу и включают в себя релаксационные процессы солей типа MgSO4, B(OH)3, MgCO3.

значение пространственного затухания \f$ \phi_{V} \f$ имеет экспотенциальную форму

\f[
\phi_{V} = exp(-\alpha s),
\f]

где \f$ s \f$ - путь пройденный лучем.

&alpha; - коэффициент затухания, который либо рассчитывается теоретически по формуле Торпа

\f[
\alpha = \frac{40 f^{2}}{4100 + f^{2}}+\frac{0.1 f^{2}}{1 + f^{2}},
\f]

f - частота в килогерцах.

Либо для данной частоты берется значение из базы.

*/