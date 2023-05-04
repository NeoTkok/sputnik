# Моделирование маневров 
_Задание для кафедры..._
### theory
#### Описание орбиты
<p>Орбита в пространстве однозначно задается пятью параметрами:<p>

__$p$__ - параметр орбиты
__$e$__ - эксцентриситет
__$\Omega$__ - долгота восходящего узла
__$\omega$__ - аргумент перигея
__$i$__ - наклонение
  
![Параметры орбиты](/images/orb.png "Параметры орбиты")

Чтобы задать положения КА в определенный момент времени - требуется еще один параметр(например средняя аномалия), однако в данной задачи не требуется точное определение этого момента, поэтому ограничемся параметром:
__$\vec {sign}$__ = $
(S_x S_y S_z)^T $ , где $\vec S = [\vec{R},\vec{V}]$.
который показывает в какую сторону движется орбита относительно её ориентации. это потребуется при изменини скорости орбиты в точках прикосновениях с целевой. Если у двух орбит $(\vec S_1,\vec S_2) > 0$ - то они "симфазны". иначе скорость совсем другую и в другом направлении придется поддавать.

Далее в кратце формулы, которые использовались для решения задачи:

1. Единичный вектор радиального направления:
$\vec {r_0}=
\left\{
 \begin{matrix}
   \cos(\Omega)\cos(\omega)-\sin(\Omega)\sin(\omega)\cos(i)\\
   \sin(\Omega)\cos(\omega)+\sin(\Omega)\sin(\omega)\cos(i)\\
   \sin(\omega)\sin(i)
  \end{matrix} 
\right\}
$
1. Уравнение Плоскости:
$Ax+By+Cz+D = 0$
1. Вектор пересечений двух плоскостей:

$\vec a = \frac{[\vec{n_1},\vec{n_2}]}{|[\vec{n_1},\vec{n_2}]|}$

4. Связи между геометрическими составляющими орбиты:
   
|        | $a$        | $p$                |$r_\pi$ |$r_a$ |
| ------ | ---------- | ------------------ | - | - |
| $a$    | $a$        |$\frac {p}{1-e^2}$|$\frac {r_\pi}{1-e}$|$\frac {r_a}{1+e}$|
| $p$    |$a(1-e^{2})$| $p$                |$r_\pi(1+e)$          |$r_a(1-e)$|
| $r_\pi$|$a(1-e)$    |$\frac {p}{1+e}$    |$r_\pi$|$r_a\frac{1-e}{1+e}$|
| $r_a$  |$a(1+e)$    |$\frac {p}{1-e}$    |$r_\pi\frac{1+e}{1-e}$|$r_a$|

5. Далее за правый полюс выберем планету, тогда справедлив:
   
   $R(\phi) = \frac{p}{1+e \cos(\phi)}$

6. Интеграл энергии:
$h = - \frac{\mu}{p}(1 -e^2)$

7. Радиальная скорость:
$V_r = \sqrt{\frac{\mu}{p}}e\sin(\phi)$

8. Трансверсальная скорость:
$V_\tau = \sqrt{\frac{\mu}{p}}(1+e\cos(\phi))$

9.  Изменение параметров орбиты при изменении скорости в перицентре на $\Delta V_\pi:$
    
    $ p' = \frac{1}{\mu}(\frac{(V_\pi+\Delta V_\pi)p}{(1+e)})^2$; $e' = \frac{(V_\pi+\Delta V_\pi)^2p}{\mu(1+e)}-1$;
     где $(p,e,V_\pi)$- стратые параметры, $\mu$ - гравитационный параметр

10. Аналогично, но только уже в апоцентре:
    
    $ p' = \frac{1}{\mu}(\frac{(V_\pi+\Delta V_\pi)p}{(1-e)})^2$; $e' =1 - \frac{(V_\pi+\Delta V_\pi)^2p}{\mu(1-e)}$;

#### Oсновные моменты:

* Библиотека состоит из 2-х _hpp_ и тестеровочного _cpp_ 
_Orbit.hpp_ : описана вспомогательная структура орбиты со множествои различных функций
_Spitnik.hpp_: описан основоной класс КА, в котором реализована главная функция _manevr()_ 
* Всё собрано в один проект с помощью _make_(автоматизированнa сборка с помощью _cmake_) 
* Основые реализационные моменты описаны в комментариях в программе
* Также имеется много тонкостей в этом плане, в плане оптимизации и логичности рассуждений, которые целесобразней обсудить лично
* Углы на эллипсе, между касательными, колличество корней в функции $f$ и $df$ при различных параметрах (см код) проверял на geogebrа (с кучей скринов и видео) - всё сошлось 

##### B кратце про manevr("орбита", "eps") - метод класса sputnik
1. Проверяет на совпадение плоскостей.
2. Если не совпали, то ищет точные точки для маневра(на прямой пересечения плоскостей орбит) - и запоминает скорость угол $\phi$ - от перицентра против ч.с. и угол под которым требуется поддать мгновенный импульс. 
3. Затем функция выводит параметры всмопогательной орбиты, на которую перешел КА
4. Далее мы уже находимся в одной плоскости целевой орбиты, где либо поддаём импульс в перецентре чтоб произошло касание орбит(с выводом ещё одной вспомогательной орбиты), либо сразу в оптимальном месте разворачиваемся под некоторой характерной скорости.