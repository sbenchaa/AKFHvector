/**
 * @fileOverview An optimized algebraic computing library to using asmjs
 * @author Sofiane Benchaa
 * @version 0.1b
 */

/* 
 * Copyright (C) 2014 Sofiane Benchaa
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


var AKFHvector = (function() {
    "use strict";

    var nbEntries = 1024;
    var size = nbEntries * 4; //1024 * Float32Array.BYTES_PER_ELEMENT, min allocation is 4096
    var buffer = new ArrayBuffer(size);
    var array = new Float64Array(buffer);

    var AKFHvectorClass = function(stdlib, userlib, heap) {
        "use asm";

        // declaration of globals
        var HEAPF64 = new stdlib.Float64Array(heap);
        var pointer = 0;

        var sqrt = stdlib.Math.sqrt;
        var abs = stdlib.Math.abs;

        /**
         * Rewinds the pointer at the beginning of heap
         *
         * @returns {undefined}
         */
        function resetHEAP() {
            pointer = 0;
        }

        /**
         * Sets a sized vector,
         * it returns a memory heap access,
         * each call shifts the heap pointer until to reach the maximum allocations
         *
         * @param {int} size
         * @returns {int}
         */
        function create(size) {
            size = size | 0;
            var lastPtr = 0;
            lastPtr = pointer;
            pointer = (pointer + size | 0) | 0; // n*4 bytes
            return lastPtr | 0;
        }

        /**
         * Gets the sign of pointed element
         *
         * @param {int} ptr
         * @param {int} index
         * @returns {float} -1.0|0.0|1.0
         */
        function sign(ptr, index) {
            ptr = ptr | 0;
            index = index | 0;
            var n = 0.0;
            n = +HEAPF64[ (ptr + index << 3) >> 3 ];
            if (n < 0.0)
                return -1.0;
            if (n > 0.0)
                return 1.0;
            return 0.0;
        }

        /**
         * Gets the value of the pointer
         *
         * @param {int} ptr
         * @param {int} index
         * @returns {LinearClass.stdlib.Float32Array}
         */
        function fetch(ptr, index) {
            ptr = ptr | 0;
            index = index | 0;
            return +HEAPF64[ (ptr + index << 3) >> 3 ];
        }

/////////// Vector2

        /**
         * Sets the xyz values of the vector
         *
         * @param {int} ptr
         * @param {float} x
         * @param {float} y
         * @returns {undefined}
         */
        function putVec2(ptr, x, y) {
            ptr = ptr | 0;
            x = +x;
            y = +y;
            HEAPF64[ (ptr + 0 << 3) >> 3 ] = x; // shift must be 2 for 32bit 3 for 64bit
            HEAPF64[ (ptr + 1 << 3) >> 3 ] = y;
            return;
        }

        function copyVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
        }

        /**
         * Clear all values of the vector,
         * filling by zero
         * x = 0.0;
         * y = 0.0;
         *
         * @param {int} ptr
         * @returns {undefined}
         */
        function resetVec2(ptr) {
            ptr = ptr | 0;
            putVec2(ptr, 0.0, 0.0);
            return;
        }

        /**
         * Tests two vectors
         * Axy == Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function equalVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axy <= Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function lequalVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axy >= Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function gequalVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0)) | 0;
        }

        /**
         * Adds two vectors
         * Axy += Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function addVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is added by a scalar
         * Axy += v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function addVec2Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + v;
            return;
        }

        /**
         * Substracts two vectors
         * Axy -= Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function subVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is substracted by a scalar
         * Axy -= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function subVec2Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - v;
            return;
        }

        /**
         * Mulitplies two vectors
         * Axy *= Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function mulVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is multiplied by a scalar
         * Axy *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function mulVec2Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * v;
            return;
        }

        /**
         * Divide two vectors
         * Axy /= Bxy
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function divVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is divided by a scalar
         * Axy *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function divVec2Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / v;
            return;
        }

        /**
         * The vector is rounded to the inferior integer,
         * beware the result is a float value
         * x = 10.5 => 10.0;
         * y = 3.14 => 3.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function floorVec2(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x);
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 - (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y);
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 - (y % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is rounded to the superior integer,
         * beware the result is a float value
         * x = 10.5 => 11.0;
         * y = 3.14 => 4.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function ceilVec2(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x) + 1.0;
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 + (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y) + 1.0;
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 + (y % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is inversed
         * x = 10.5 => -10.5;
         * y = 3.14 => -3.14;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function negateVec2(ptr1) {
            ptr1 = ptr1 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 0 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 1 << 3) >> 3 ]);
            return;
        }

        /**
         * Realizes a dot product between two vectors
         * Axy . Bxy
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function dotVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var ab0 = 0.0;
            var ab1 = 0.0;
            var ab2 = 0.0;
            ab0 = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            ab1 = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            return +(ab0 + ab1);
        }

        /**
         * Realizes a dot product between two vectors
         * Axy ^ Bxy
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function crossVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var ax = 0.0;
            var ay = 0.0;
            var bx = 0.0;
            var by = 0.0;
            ax = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            ay = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            bx = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            by = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ];

            return +(ax * by - ay * bx);
        }

        /**
         * Interpolates two vectors by a scalar
         * Axy += (Bxy - Axy) * v
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @param {type} v
         * @returns {undefined}
         */
        function lerpVec2(ptr1, ptr2, v) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            v = +v;
            var bSa0 = 0.0;
            var bSa1 = 0.0;
            bSa0 = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            bSa1 = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + bSa0 * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + bSa1 * v;
            return;
        }

        /**
         * Computes the magnitude of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthSqVec2(ptr) {
            ptr = ptr | 0;
            return +dotVec2(ptr, ptr);
        }

        /**
         * Computes the length of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthVec2(ptr) {
            ptr = ptr | 0;
            return +sqrt(+dotVec2(ptr, ptr));
        }

        /**
         * Computes the Euclidian distance between two vectors
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function distanceVec2(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            return +lengthVec2(ptr1);
        }

        /**
         * Computes the Manhattan distance between two vectors
         *
         * @param {type} ptr
         * @returns {unresolved}
         */
        function lengthManVec2(ptr) {
            ptr = ptr | 0;
            return +(+abs(+HEAPF64[ (ptr + 0 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 1 << 3) >> 3 ]));
        }

        /**
         * Normalizes the vector,
         * reduces the vector between [0.0, 1.0]
         *
         * @param {type} ptr
         * @returns {undefined}
         */
        function normalizeVec2(ptr) {
            ptr = ptr | 0;
            divVec2Scalar(ptr, +lengthVec2(ptr));
            return;
        }

        /////////// Vector3

        /**
         * Sets the xyz values of the vector
         *
         * @param {int} ptr
         * @param {float} x
         * @param {float} y
         * @param {float} z
         * @returns {undefined}
         */
        function putVec3(ptr, x, y, z) {
            ptr = ptr | 0;
            x = +x;
            y = +y;
            z = +z;
            HEAPF64[ (ptr + 0 << 3) >> 3 ] = x; // shift must be 2 for 32bit 3 for 64bit
            HEAPF64[ (ptr + 1 << 3) >> 3 ] = y;
            HEAPF64[ (ptr + 2 << 3) >> 3 ] = z;
            return;
        }

        function copyVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
        }

        /**
         * Clear all values of the vector,
         * filling by zero
         * x = 0.0;
         * y = 0.0;
         * z = 0.0;
         *
         * @param {int} ptr
         * @returns {undefined}
         */
        function resetVec3(ptr) {
            ptr = ptr | 0;
            putVec3(ptr, 0.0, 0.0, 0.0);
            return;
        }

        /**
         * Tests two vectors
         * Axyz == Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function equalVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axyz <= Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function lequalVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axyz >= Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function gequalVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0)) | 0;
        }

        /**
         * Adds two vectors
         * Axyz += Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function addVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is added by a scalar
         * Axyz += v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function addVec3Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + v;
            return;
        }

        /**
         * Substracts two vectors
         * Axyz -= Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function subVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] - HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is substracted by a scalar
         * Axyz -= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function subVec3Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] - v;
            return;
        }

        /**
         * Mulitplies two vectors
         * Axyz *= Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function mulVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is multiplied by a scalar
         * Axyz *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function mulVec3Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * v;
            return;
        }

        /**
         * Divide two vectors
         * Axyz /= Bxyz
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function divVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] / HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is divided by a scalar
         * Axyz *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function divVec3Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] / v;
            return;
        }

        /**
         * The vector is rounded to the inferior integer,
         * beware the result is a float value
         * x = 10.5 => 10.0;
         * y = 3.14 => 3.0;
         * z = 0.9 => 0.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function floorVec3(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            var z = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            z = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x);
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 - (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y);
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 - (y % 1.0))) - 1.0));
            if (z > 0.0)
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~z);
            else
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~(z - +(~~(-1.0 - (z % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is rounded to the superior integer,
         * beware the result is a float value
         * x = 10.5 => 11.0;
         * y = 3.14 => 4.0;
         * z = 0.9 => 1.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function ceilVec3(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            var z = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            z = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x) + 1.0;
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 + (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y) + 1.0;
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 + (y % 1.0))) - 1.0));
            if (z > 0.0)
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~z) + 1.0;
            else
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~(z - +(~~(-1.0 + (z % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is inversed
         * x = 10.5 => -10.5;
         * y = 3.14 => -3.14;
         * z = 0.9 => -0.9;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function negateVec3(ptr1) {
            ptr1 = ptr1 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 0 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 1 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 2 << 3) >> 3 ]);
            return;
        }

        /**
         * Realizes a dot product between two vectors
         * Axyz . B.xyz
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function dotVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var ab0 = 0.0;
            var ab1 = 0.0;
            var ab2 = 0.0;
            ab0 = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            ab1 = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            ab2 = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            return +(ab0 + ab1 + ab2);
        }

        /**
         * Realizes a dot product between two vectors
         * Axyz ^ B.xyz
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {undefined}
         */
        function crossVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var ax = 0.0;
            var ay = 0.0;
            var az = 0.0;
            var bx = 0.0;
            var by = 0.0;
            var bz = 0.0;
            ax = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            ay = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            az = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            bx = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            by = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            bz = +HEAPF64[ (ptr2 + 2 << 3) >> 3 ];

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(ay * bz - az * by);
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(az * bx - ax * bz);
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(ax * by - ay * bx);
            return;
        }

        /**
         * Interpolates two vectors by a scalar
         * Axyz += (Bxyz - Axyz) * v
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @param {type} v
         * @returns {undefined}
         */
        function lerpVec3(ptr1, ptr2, v) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            v = +v;
            var bSa0 = 0.0;
            var bSa1 = 0.0;
            var bSa2 = 0.0;
            bSa0 = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            bSa1 = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            bSa2 = +HEAPF64[ (ptr2 + 2 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + bSa0 * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + bSa1 * v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + bSa2 * v;
            return;
        }

        /**
         * Computes the magnitude of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthSqVec3(ptr) {
            ptr = ptr | 0;
            return +dotVec3(ptr, ptr);
        }

        /**
         * Computes the length of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthVec3(ptr) {
            ptr = ptr | 0;
            return +sqrt(+dotVec3(ptr, ptr));
        }

        /**
         * Computes the Euclidian distance between two vectors
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function distanceVec3(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 2 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            return +lengthVec3(ptr1);
        }

        /**
         * Computes the Manhattan distance between two vectors
         *
         * @param {type} ptr
         * @returns {unresolved}
         */
        function lengthManVec3(ptr) {
            ptr = ptr | 0;
            return +(+abs(+HEAPF64[ (ptr + 0 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 1 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 2 << 3) >> 3 ]));
        }

        /**
         * Normalizes the vector,
         * reduces the vector between [0.0, 1.0]
         *
         * @param {type} ptr
         * @returns {undefined}
         */
        function normalizeVec3(ptr) {
            ptr = ptr | 0;
            divVec3Scalar(ptr, +lengthVec3(ptr));
            return;
        }

        /////////// Vector4

        /**
         * Sets the xyzw values of the vector
         *
         * @param {int} ptr
         * @param {float} x
         * @param {float} y
         * @param {float} z
         * @param {float} w
         * @returns {undefined}
         */
        function putVec4(ptr, x, y, z, w) {
            ptr = ptr | 0;
            x = +x;
            y = +y;
            z = +z;
            w = +w;
            HEAPF64[ (ptr + 0 << 3) >> 3 ] = x; // shift must be 2 for 32bit 3 for 64bit
            HEAPF64[ (ptr + 1 << 3) >> 3 ] = y;
            HEAPF64[ (ptr + 2 << 3) >> 3 ] = z;
            HEAPF64[ (ptr + 3 << 3) >> 3 ] = w;
            return;
        }

        function copyVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
        }

        /**
         * Clear all values of the vector,
         * filling by zero
         * x = 0.0;
         * y = 0.0;
         * z = 0.0;
         * w = 0.0;
         *
         * @param {int} ptr
         * @returns {undefined}
         */
        function resetVec4(ptr) {
            ptr = ptr | 0;
            putVec4(ptr, 0.0, 0.0, 0.0, 0.0);
            return;
        }

        /**
         * Tests two vectors
         * Axyzw == Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function equalVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            var d = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            d = (+HEAPF64[ (ptr1 + 3 << 3) >> 3 ] == +HEAPF64[ (ptr2 + 3 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0) & (d | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axyzw <= Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function lequalVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            var d = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            d = (+HEAPF64[ (ptr1 + 3 << 3) >> 3 ] <= +HEAPF64[ (ptr2 + 3 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0) & (d | 0)) | 0;
        }

        /**
         * Tests two vectors
         * Axyzw >= Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {int}
         */
        function gequalVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var a = 0;
            var b = 0;
            var c = 0;
            var d = 0;
            a = (+HEAPF64[ (ptr1 + 0 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 0 << 3) >> 3 ]) | 0;
            b = (+HEAPF64[ (ptr1 + 1 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 1 << 3) >> 3 ]) | 0;
            c = (+HEAPF64[ (ptr1 + 2 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 2 << 3) >> 3 ]) | 0;
            d = (+HEAPF64[ (ptr1 + 3 << 3) >> 3 ] >= +HEAPF64[ (ptr2 + 3 << 3) >> 3 ]) | 0;
            return  ((a | 0) & (b | 0) & (c | 0) & (d | 0)) | 0;
        }

        /**
         * Adds two vectors
         * Axyzw += Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function addVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = HEAPF64[ (ptr1 + 3 << 3) >> 3 ] + HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is added by a scalar
         * Axyzw += v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function addVec4Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + v;
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] + v;
            return;
        }

        /**
         * Substracts two vectors
         * Axyzw -= Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function subVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] - HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = HEAPF64[ (ptr1 + 3 << 3) >> 3 ] - HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is substracted by a scalar
         * Axyzw -= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function subVec4Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] - v;
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] - v;
            return;
        }

        /**
         * Mulitplies two vectors
         * Axyzw *= Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function mulVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = HEAPF64[ (ptr1 + 3 << 3) >> 3 ] * HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is multiplied by a scalar
         * Axyzw *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function mulVec4Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * v;
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] * v;
            return;
        }

        /**
         * Divide two vectors
         * Axyzw /= Bxyzw
         *
         * @param {int} ptr1
         * @param {int} ptr2
         * @returns {undefined}
         */
        function divVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = HEAPF64[ (ptr1 + 2 << 3) >> 3 ] / HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = HEAPF64[ (ptr1 + 3 << 3) >> 3 ] / HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
            return;
        }

        /**
         * The vector is divided by a scalar
         * Axyzw *= v
         *
         * @param {int} ptr1
         * @param {float} v
         * @returns {undefined}
         */
        function divVec4Scalar(ptr1, v) {
            ptr1 = ptr1 | 0;
            v = +v;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] / v;
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] / v;
            return;
        }

        /**
         * The vector is rounded to the inferior integer,
         * beware the result is a float value
         * x = 10.5 => 10.0;
         * y = 3.14 => 3.0;
         * z = 0.9 => 0.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function floorVec4(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            var z = 0.0;
            var w = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            z = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            w = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x);
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 - (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y);
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 - (y % 1.0))) - 1.0));
            if (z > 0.0)
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~z);
            else
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~(z - +(~~(-1.0 - (z % 1.0))) - 1.0));
            if (w > 0.0)
                HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +(~~w);
            else
                HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +(~~(w - +(~~(-1.0 - (w % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is rounded to the superior integer,
         * beware the result is a float value
         * x = 10.5 => 11.0;
         * y = 3.14 => 4.0;
         * z = 0.9 => 1.0;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function ceilVec4(ptr1) {
            ptr1 = ptr1 | 0;
            var x = 0.0;
            var y = 0.0;
            var z = 0.0;
            var w = 0.0;
            x = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            y = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            z = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            w = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ];
            if (x > 0.0)
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~x) + 1.0;
            else
                HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +(~~(x - +(~~(-1.0 + (x % 1.0))) - 1.0));
            if (y > 0.0)
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~y) + 1.0;
            else
                HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +(~~(y - +(~~(-1.0 + (y % 1.0))) - 1.0));
            if (z > 0.0)
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~z) + 1.0;
            else
                HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +(~~(z - +(~~(-1.0 + (z % 1.0))) - 1.0));
            if (w > 0.0)
                HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +(~~w) + 1.0;
            else
                HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +(~~(w - +(~~(-1.0 + (w % 1.0))) - 1.0));
            return;
        }

        /**
         * The vector is inversed
         * x = 10.5 => -10.5;
         * y = 3.14 => -3.14;
         * z = 0.9 => -0.9;
         *
         * @param {type} ptr1
         * @returns {undefined}
         */
        function negateVec4(ptr1) {
            ptr1 = ptr1 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 0 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 1 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 2 << 3) >> 3 ]);
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = -(+HEAPF64[ (ptr1 + 3 << 3) >> 3 ]);
            return;
        }

        /**
         * Realizes a dot product between two vectors
         * Axyzw . B.xyzw
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function dotVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            var ab0 = 0.0;
            var ab1 = 0.0;
            var ab2 = 0.0;
            var ab3 = 0.0;
            ab0 = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 0 << 3) >> 3 ];
            ab1 = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 1 << 3) >> 3 ];
            ab2 = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 2 << 3) >> 3 ];
            ab3 = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] * +HEAPF64[ (ptr2 + 3 << 3) >> 3 ];
            return +(ab0 + ab1 + ab2 + ab3);
        }

        /**
         * Interpolates two vectors by a scalar
         * Axyzw += (Bxyzw - Axyzw) * v
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @param {type} v
         * @returns {undefined}
         */
        function lerpVec4(ptr1, ptr2, v) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            v = +v;
            var bSa0 = 0.0;
            var bSa1 = 0.0;
            var bSa2 = 0.0;
            var bSa3 = 0.0;
            bSa0 = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            bSa1 = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            bSa2 = +HEAPF64[ (ptr2 + 2 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            bSa3 = +HEAPF64[ (ptr2 + 3 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 3 << 3) >> 3 ];

            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 0 << 3) >> 3 ] + bSa0 * v;
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 1 << 3) >> 3 ] + bSa1 * v;
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 2 << 3) >> 3 ] + bSa2 * v;
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr1 + 3 << 3) >> 3 ] + bSa3 * v;
            return;
        }

        /**
         * Computes the magnitude of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthSqVec4(ptr) {
            ptr = ptr | 0;
            return +dotVec4(ptr, ptr);
        }

        /**
         * Computes the length of a vector
         *
         * @param {type} ptr
         * @returns {float}
         */
        function lengthVec4(ptr) {
            ptr = ptr | 0;
            return +sqrt(+dotVec4(ptr, ptr));
        }

        /**
         * Computes the Euclidian distance between two vectors
         *
         * @param {type} ptr1
         * @param {type} ptr2
         * @returns {float}
         */
        function distanceVec4(ptr1, ptr2) {
            ptr1 = ptr1 | 0;
            ptr2 = ptr2 | 0;
            HEAPF64[ (ptr1 + 0 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 0 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 0 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 1 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 1 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 1 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 2 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 2 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 2 << 3) >> 3 ];
            HEAPF64[ (ptr1 + 3 << 3) >> 3 ] = +HEAPF64[ (ptr2 + 3 << 3) >> 3 ] - +HEAPF64[ (ptr1 + 3 << 3) >> 3 ];
            return +lengthVec4(ptr1);
        }

        /**
         * Computes the Manhattan distance between two vectors
         *
         * @param {type} ptr
         * @returns {unresolved}
         */
        function lengthManVec4(ptr) {
            ptr = ptr | 0;
            return +(+abs(+HEAPF64[ (ptr + 0 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 1 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 2 << 3) >> 3 ]) +
                    +abs(+HEAPF64[ (ptr + 3 << 3) >> 3 ]));
        }

        /**
         * Normalizes the vector,
         * reduces the vector between [0.0, 1.0]
         *
         * @param {type} ptr
         * @returns {undefined}
         */
        function normalizeVec4(ptr) {
            ptr = ptr | 0;
            divVec4Scalar(ptr, +lengthVec4(ptr));
            return;
        }

        return {
            resetHEAP: resetHEAP,
            create: create,
            fetch: fetch,
            sign: sign,
            putVec2: putVec2,
            copyVec2: copyVec2,
            resetVec2: resetVec2,
            equalVec2: equalVec2,
            lequalVec2: lequalVec2,
            gequalVec2: gequalVec2,
            addVec2: addVec2,
            addVec2Scalar: addVec2Scalar,
            subVec2: subVec2,
            subVec2Scalar: subVec2Scalar,
            mulVec2: mulVec2,
            mulVec2Scalar: mulVec2Scalar,
            divVec2: divVec2,
            divVec2Scalar: divVec2Scalar,
            floorVec2: floorVec2,
            ceilVec2: ceilVec2,
            negateVec2: negateVec2,
            dotVec2: dotVec2,
            lerpVec2: lerpVec2,
            crossVec2: crossVec2,
            lengthSqVec2: lengthSqVec2,
            lengthVec2: lengthVec2,
            lengthManVec2: lengthManVec2,
            distanceVec2: distanceVec2,
            normalizeVec2: normalizeVec2,
            putVec3: putVec3,
            copyVec3: copyVec3,
            resetVec3: resetVec3,
            equalVec3: equalVec3,
            lequalVec3: lequalVec3,
            gequalVec3: gequalVec3,
            addVec3: addVec3,
            addVec3Scalar: addVec3Scalar,
            subVec3: subVec3,
            subVec3Scalar: subVec3Scalar,
            mulVec3: mulVec3,
            mulVec3Scalar: mulVec3Scalar,
            divVec3: divVec3,
            divVec3Scalar: divVec3Scalar,
            floorVec3: floorVec3,
            ceilVec3: ceilVec3,
            negateVec3: negateVec3,
            dotVec3: dotVec3,
            lerpVec3: lerpVec3,
            crossVec3: crossVec3,
            lengthSqVec3: lengthSqVec3,
            lengthVec3: lengthVec3,
            lengthManVec3: lengthManVec3,
            distanceVec3: distanceVec3,
            normalizeVec3: normalizeVec3,
            putVec4: putVec4,
            copyVec4: copyVec4,
            resetVec4: resetVec4,
            equalVec4: equalVec4,
            lequalVec4: lequalVec4,
            gequalVec4: gequalVec4,
            addVec4: addVec4,
            addVec4Scalar: addVec4Scalar,
            subVec4: subVec4,
            subVec4Scalar: subVec4Scalar,
            mulVec4: mulVec4,
            mulVec4Scalar: mulVec4Scalar,
            divVec4: divVec4,
            divVec4Scalar: divVec4Scalar,
            floorVec4: floorVec4,
            ceilVec4: ceilVec4,
            negateVec4: negateVec4,
            dotVec4: dotVec4,
            lerpVec4: lerpVec4,
            lengthSqVec4: lengthSqVec4,
            lengthVec4: lengthVec4,
            lengthManVec4: lengthManVec4,
            distanceVec4: distanceVec4,
            normalizeVec4: normalizeVec4
        };
    };


    var akfh = AKFHvectorClass(window, {}, buffer);
    akfh.count = 0;
    /**
     * Gets an array of vector
     *
     * @param {int} ptr
     * @param {int} len
     * @returns {Float32Array}
     */
    akfh.getVector = function(ptr, len) {
        return array.subarray(ptr, ptr + len);
    };

    /**
     *
     * @param {type} size
     * @returns {undefined}
     */
    akfh.clearHEAP = function(size) {
        akfh.count = 0;
        akfh.resetHEAP( );

    };

    /**
     *
     * @param {type} size
     * @returns {Number|@exp;_L5@pro;nbEntries|@exp;_L5@pro;size|@exp;_L5@pro;AKFHvectorClass@pro;pointer|int}
     */
    akfh.allocate = function(size) {
        if (akfh.count + size > nbEntries) {
            console.error("The maximum of entries was reached ", akfh.count + size, ". Limit is ", nbEntries);
        } else {
            akfh.count += size;
            return akfh.create(size);
        }
    };

    return akfh;
})();
