AKFHvector*
==========

An optimized algebraic computing library using asmjs.


EN: The TypedArray are essential in the buffers managing. Because they are typed,
they play an important role in the research off efficient unlike to the simple arrays.
Eg: the webGL buffers are required for building and transfering datas.

This library allow usual algebraic operations like the vectors, matrices and many things other.

In the practice, it looks like a glMatrix, Sylvester or other javascript math library but not only, 
indeed it has an optimal approach. 
It embeds advanced features in order to provide the best API through ASMJS and other in the future.

--

FR: Les TypedArray sont les éléments incontournables dans la gestion des opérations
en buffer. Parce qu'ils sont typés, ils permettent un gain en performance significatif
par rapport au simple tableau javascript. Ils sont par exemple requis en WebGL pour ses buffers.

Cette bibliothèque offre les outils nécessaires aux calcules mathématiques. 
Il est comparable en ce point à glMatrix, Sylvester, etc. mais pas seulement, car il
dispose d'une approche très optimisée. En effet, il tire pleinement partie des avantages
techniques des navigateurs. Il est dans sa première mouture écrit en ASMJS.

AKFHvector affiche ainsi de très hautes performances.

*Al Kafi-Fil-Hisab (which is sufficient for the calculation)

--

```javascript
//creates and instantiates a sizing vector 3. An uint value is attributed to myVec3 like pointer.
var myVec3 = AKFHvector32.allocate(3);

//myVec3 is setted by this following values [0.6, 0.7, 0.8].
AKFHvector32.putVec3(myVec3, 0.6, 0.7, 0.8);

//An addition on myVec3 is performed by this following values [0.4, 1.0, 2.2]. It equals now at [1.0, 1.7, 3.0].
AKFHvector32.addVec3(myVec3, 0.4, 1.0, 2.2);
```
