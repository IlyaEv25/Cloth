/**
    Главный класс симуляции. Все взаимодействия с параметрами ткани через него.
    
 */
/** Главный класс симуляции. Все взаимодействия с параметрами ткани через него.  
         * @property {Integer} size Нечетное число! Количество вершин в одном ряду и количество рядов.
         * @property {WebGLContext} gl Указатель на контекст WebGL. 
         * @property {Float[]} X Массив для буфера вершин. В нем содержатся координаты вершин и информация для наложения текстур.
         * @property {Float[]} Y Массив положений вершин в предыдуший момент времени.
         * @property {Vector2D[]} u Массив 2-мерных векторов, обозначающих начальное положение вершин.
         * @property {Vector3D[]} w_u Массив 3-мерных векторов, обозначающих производную линейного отображения по направлению u для одного треугольника.(уравнение 1)
         * @property {Vector3D[]} w_v Массив 3-мерных векторов, обозначающих производную линейного отображения по направлению v для одного треугольника. (уравнение 1)
         * @property {Vector3D[]} F Массив 3-мерных векторов, обозначающих силы, действующие на вершины.
         * @property {Vector3D[]} V Массив 3-мерных векторов, обозначающих скорости вершин.
         * @property {Vector3D[]} z Массив 3-мерных векторов, обозначающих ограничения на изменение скорости вершин.
         * @property {Matrix2D[]} Uinv Массив 2-мерных матриц, по которым из координат в конкретный момент строятся производные линейных отображений.(уравнение 3)
         * @property {Matrix3D[]} S Массив 3-мерных матриц, которые реализуют ограничения движений вершин.
         * 
         * @property {Vector3D[]} triangles Массив 3-мерных векторов, в каждом из которых - индексы трех вершин, входящих в конкретный треугольник.
         * @property {Integer[]} tri_inds Массив индексов вершин. Плоская версия массива triangles, нужна для использования в буффере индексов.
         * @property {Vector2D[]} tri_pairs Массив пар индексов треугольников.
         * 
         * @property {Float} minv - Обратная масса одной вершины.
         * @property {Float} a В данной реализации 1. Изначально предполагалось, что это велиична, равная площади треугольника в начальном положении. 
         * @property {Float} g Гравитационная постоянная.
         * @property {Float} ks Коэффициент перед величинами, описывающими растяжение.
         * @property {Float} kd Коэффициент перед величинами, описывающими дэмпинг сил, связанных с растяжением.
         * @property {Float} kb Коэффициент перед величинами, описывающими изгиб.
         * @property {Float} kdb Коэффициент перед величинами, описывающими дэмпинг сил, связанных с изгибом.
         * @property {Float} bu Величина, к которой в динамике стремится длина производной по направлению u(процедура минимизирует функцию (|w_u| - bu)).
         * @property {Float} bv Величина, к которой в динамике стремится длина производной по направлению v(процедура минимизирует функцию (|w_v| - bv)).
         * 
         * @property {BufferObject} vbo Буфер вершин.
         * @property {BuffetObject} ebo Буфер индексов.
         */

class Cloth
{
    /**
     * 
     * @param {Integer} size Нечетное число! Количество вершин в одном ряду и количество рядов.
     * @param {WebGLcontext} gl Указатель на контекст WebGL. 
     */
    
    constructor(size, gl)
    {

        this.size = size;
        this.gl = gl;
        this.X = new Float32Array(size * size * 5);
        this.Y = new Float32Array(size * size * 3);
        this.u = [];
        this.w_u = [];
        this.w_v = [];
        this.F = [];
        this.V = [];
        this.z = [];
        this.Uinv = [];
        this.S = [];


        this.triangles = [];
        this.tri_inds = [];
        this.tri_pairs = [];


        this.minv = 1;
        this.a = 1.0;
        this.g = 9.8;
        this.ks = 1000.0;
        this.kd = 10.8;
        this.kb = 40.0;
        this.kdb = 18.8;
        this.bu = 1;
        this.bv = 1;

        let scaling = 20; // Умножает все координаты на соответствующее число. Делает модель стабильнее.

        /**
         * Триангуляция поверхности ткани. Выглядит следующим образом:                          |/|\|/|\|
         * Одновременно инцииализируем начальное положение, положение с течением времени,       |\|/|\|/|
         * координаты для дальнейшего наложения текстур, начальные силы, параметры(z)
         * и матрицу(S), которые будут использоваться для наложения ограничений на изменение скорости(dv).
         */

        for (let i = 0; i < size; i++)
        {
            for (let j = 0; j < size; j++)
            {
                this.u.push([scaling * j/(size - 1), scaling * i/(size - 1)]);
                this.X[5*(i * size + j)] = scaling * j/(size - 1);
                this.X[5*(i * size + j) + 1] = scaling * i/(size - 1);
                this.X[5*(i * size + j) + 2] = 0;
                this.X[5*(i * size + j) + 3] = j/(size - 1);
                this.X[5*(i * size + j) + 4] = i/(size - 1);

                this.Y[3*(i * size + j)] = scaling * j/(size - 1);
                this.Y[3*(i * size + j) + 1] = scaling * i/(size - 1);
                this.Y[3*(i * size + j) + 2] = 0;

                this.F.push([0, 0, 0]);
                this.V.push([0, 0, 0]);
                this.z.push([0, 0, 0]);
                
                if (i == 0)
                    this.S.push([[0, 0, 0], [0, 0, 0], [0, 0, 0]]); 
                else
                    this.S.push([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
            }                

        }
        
        /**
         * Инициализируем массив, каждый элемент которого - индексы трех точек в большом массиве Х,
         * которые образуют треугольник.
         */
            
        for (let j = 0; j < size - 1; j++)
        {
            if (j % 2 == 0)
            {
                for (let k = j * size; k < j*size + size - 2; k = k + 2)
                {
                    this.triangles.push([k, k + 1, k + size]);
                    this.triangles.push([k + 1, k + size, k + size + 1]);
                    this.triangles.push([k + 1, k + size + 1, k + size + 2]);
                    this.triangles.push([k + 1, k + 2, k + size + 2]);
                    
                }
            }
            else
            {
                for (let k = j * size; k < j*size + size - 2; k = k + 2)
                {
                    this.triangles.push([k, k + size, k + size + 1]);
                    this.triangles.push([k, k + 1, k + size + 1]);
                    this.triangles.push([k + 1, k + 2, k + size + 1]);
                    this.triangles.push([k + 2, k + size + 1, k + size + 2]);
                    
                }
            }        
        }

        /**
         * Инициализируем массив индексов вершин тругольников, массив матриц обратных U,
         * массив производных линейных отображений, характеризующих растяжение треугольников,
         * массив пар треугольников.
         */

        let NTR = 2 * (size - 1); // количество треугольников в ряде
        for (let k = 0; k < this.triangles.length; k++)
        {
            this.tri_inds.push(this.triangles[k][0]);
            this.tri_inds.push(this.triangles[k][1]);
            this.tri_inds.push(this.triangles[k][2]);
            this.Uinv.push(this.Umat(k));
            this.w_u.push([1, 0, 0]);
            this.w_v.push([0, 1, 0]);

            if (k % NTR != NTR - 1)
                this.tri_pairs.push([k, k + 1]);
            if ((k / NTR >> 0) % 2 == 0)
            {
                if (((k % NTR) % 4 == 1) && ((k / NTR >> 0)!= size - 2))
                {
                    this.tri_pairs.push([k, k + NTR]);
                    this.tri_pairs.push([k + 1, k + 1 + NTR]);
                }
            }
            else
            {
                if (((k % NTR) % 4 == 0) && ((k / NTR >> 0)!= size - 2))
                {
                    this.tri_pairs.push([k, k + NTR]);
                    this.tri_pairs.push([k + 3, k + 3 + NTR]);
                }
            }
        } 


        this.vbo = gl.createBuffer();
        this.ebo = gl.createBuffer();

        /**
         * Рутина WebGL. Инициализируем вершинные буфферные объекты, индексные,
         * и определяем аттрибуты вершинного массива.
         */

        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo); 
        gl.bufferData(gl.ARRAY_BUFFER, this.X, gl.STREAM_DRAW);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.ebo);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.tri_inds), gl.STATIC_DRAW);

        gl.enableVertexAttribArray(0);
        gl.enableVertexAttribArray(1);
        gl.vertexAttribPointer(0, 3, gl.FLOAT, false, this.X.BYTES_PER_ELEMENT * 5, 0);
        gl.vertexAttribPointer(1, 2, gl.FLOAT, false, this.X.BYTES_PER_ELEMENT * 5, this.X.BYTES_PER_ELEMENT * 3);

    }

    /**
     * Обнуляем силы. Перед рассчетом сил в следующий момент времени должны их обнулить.
     */
    zero()
    {
        let n;
        
        n = this.F.length;
        for (let k = 0; k < n; k++)
            this.F[k] = [0, 0, 0];

    }

    /**
     * Вычисляем матрицу U.
     * 
     * @param {Integer} triangle Индекс треугольника.
     * @returns {Matrix2D}
     */

    Umat(triangle)
    {
        let i, j, k;
        let det, du1, du2, dv1, dv2;
        i = this.triangles[triangle][0];
        j = this.triangles[triangle][1];
        k = this.triangles[triangle][2];

        du1 = this.u[j][0] - this.u[i][0];
        du2 = this.u[k][0] - this.u[i][0];
        dv1 = this.u[j][1] - this.u[i][1];
        dv2 = this.u[k][1] - this.u[i][1];

        det = du1 * dv2 - dv1 * du2;

        return [dv2/det, -du2/det, -dv1/det, du1/det];
    }

    /**
     * Считаем силы, ассоциированные с энергией, зависящей от пары треугольников.
     * Т.е. силу, сопротивляющуюся изгибу.
     * @param {Integer} pair 
     */
    calcTriPair(pair)
    {
        let i, j, k, m, l;
        let inters = [];
        let ex = [];

        let kb = this.kb;
        let kdb = this.kdb;

        /**
         * Выбираем 4 неповторяющиеся вершины пары треугольников с общим ребром таким образом,
         * чтобы j, k были вершинами на общем ребре, а i, m - вершинами, которые принадлежат только 
         * одному треугольнику в паре. 
         */
        for (let t = 0; t < 3; t++)
        {
            l = 0;
            for(let r = 0; r < 3; r++)
            {
                if (this.triangles[this.tri_pairs[pair][0]][t] == this.triangles[this.tri_pairs[pair][1]][r])
                {
                    inters.push(this.triangles[this.tri_pairs[pair][0]][t]);
                    l++;
                }
            }
            if (l == 0)
                ex.push(this.triangles[this.tri_pairs[pair][0]][t]);
        }
        
        for (let t = 0; t < 3; t++)
        {
            l = 0;
            for(let r = 0; r < 2; r++)
            {
                if (this.triangles[this.tri_pairs[pair][1]][t] == inters[r])
                    l++;
            }
            if (l == 0)
                ex.push(this.triangles[this.tri_pairs[pair][1]][t]);
        }

        i = ex[0];
        j = inters[0];
        k = inters[1];
        m = ex[1];

        /**
         * Вычисления, описанные в уравнениях 9, 12. После вычисляется вклад данной пары в силы,
         * действующие на вершины i, j, k, m.
         */
        let dx1 = [this.X[j*5] - this.X[i*5], this.X[j*5 + 1] - this.X[i*5 + 1], this.X[j*5 + 2] - this.X[i*5 + 2]];
        let dx2 = [this.X[k*5] - this.X[i*5], this.X[k*5 + 1] - this.X[i*5 + 1], this.X[k*5 + 2] - this.X[i*5 + 2]];
        let dx3 = [this.X[k*5] - this.X[j*5], this.X[k*5 + 1] - this.X[j*5 + 1], this.X[k*5 + 2] - this.X[j*5 + 2]];
        let dx4 = [this.X[m*5] - this.X[j*5], this.X[m*5 + 1] - this.X[j*5 + 1], this.X[m*5 + 2] - this.X[j*5 + 2]];

        let N1 = cross(dx1, dx2);
        let N2 = cross(dx3, dx4);
        if (dot(N1, N2) < 0)
            N2 = numxvec(-1, N2);

        let Cb = 1 - (1.0/length(N1)) * (1.0/length(N2)) * dot(N1, N2);           

        /**
         * Предполагая, что нормали почти постоянны по длине(что соответствует происходящему в симуляции),
         * пренебрегаем производными от их длины. Матрицы, соответствующие производным по нормалям, 
         * являются матрицами, ассоциированными с некоторыми векторами в смысле векторного произведения, т.е. если 
         * такую матрицу умножить на другой вектор, то это равносильно векторному произведению ассоциированного с матрицей вектора
         * с вектором, на который производится умножение.
         */
        let DCbDxi = numxvec((-1.0/length(N1)) * (1.0/length(N2)), cross(difvecs(dx2,dx1), N2));
        let DCbDxj = numxvec((-1.0/length(N1)) * (1.0/length(N2)), sumvecs(cross(numxvec(-1, dx2), N2), cross(numxvec(-1, difvecs(dx3, dx4)), N1)));
        let DCbDxk = numxvec((-1.0/length(N1)) * (1.0/length(N2)), sumvecs(cross(dx1, N2), cross(numxvec(-1, dx4), N1)));
        let DCbDxm = numxvec((-1.0/length(N1)) * (1.0/length(N2)), cross(dx3, N1));


        let dotCb = dot(DCbDxi, this.V[i]) + dot(DCbDxj, this.V[j]) + dot(DCbDxk, this.V[k]) + dot(DCbDxm, this.V[m]);

        
        let fi = numxvec(-kb * Cb - kdb * dotCb, DCbDxi); 
        let fj = numxvec(-kb * Cb - kdb * dotCb, DCbDxj);
        let fk = numxvec(-kb * Cb - kdb * dotCb, DCbDxk);
        let fm = numxvec(-kb * Cb - kdb * dotCb, DCbDxm);
        
        this.F[i] = sumvecs(this.F[i], fi);
        this.F[j] = sumvecs(this.F[j], fj);
        this.F[k] = sumvecs(this.F[k], fk);
        this.F[m] = sumvecs(this.F[m], fm);

    }

    /**
     * Считаем силы, ассоциированные с энергией, зависящей от одного треугольника.
     * Т.е. силы, сопротивляющиеся растяжению и сдвигу.
     * @param {Integer} triangle 
     */

    calcTriangle(triangle)
    {
        let i = this.triangles[triangle][0];
        let j = this.triangles[triangle][1];
        let k = this.triangles[triangle][2];

        let ks = this.ks;
        let kd = this.kd;

        /**
         * Здесь производятся вычисления, описанные в уравнениях 7, 8, 10, 11.
         */

        let dx1 = [this.X[5*j] - this.X[5*i], this.X[5*j + 1] - this.X[5*i + 1], this.X[5*j + 2] - this.X[5*i + 2]];
        let dx2 = [this.X[5*k] - this.X[5*i], this.X[5*k + 1] - this.X[5*i + 1], this.X[5*k + 2] - this.X[5*i + 2]];

        let wu = sumvecs(numxvec(this.Uinv[triangle][0], dx1), numxvec(this.Uinv[triangle][2], dx2));
        let wv = sumvecs(numxvec(this.Uinv[triangle][1], dx1), numxvec(this.Uinv[triangle][3], dx2));
        let wul = length(wu);
        let wvl = length(wv);

        let Cstru = this.a * (wul - this.bu);
        let Cstrv = this.a * (wvl - this.bv);
        let Csh =  this.a * dot(wu, wv); 

        let DwuDxi = - this.Uinv[triangle][0] - this.Uinv[triangle][2]; // Это диагональные матрицы, операции с которыми сводятся к операциям с числами на диагонали.
        let DwuDxj = this.Uinv[triangle][0];
        let DwuDxk = this.Uinv[triangle][2];
        let DwvDxi = - this.Uinv[triangle][1] - this.Uinv[triangle][3];
        let DwvDxj = this.Uinv[triangle][1];
        let DwvDxk = this.Uinv[triangle][3];

        let DCstruDxi = numxvec(this.a * 1.0/wul * DwuDxi, wu);
        let DCstruDxj = numxvec(this.a * 1.0/wul * DwuDxj, wu);
        let DCstruDxk = numxvec(this.a * 1.0/wul * DwuDxk, wu);

        let DCstrvDxi = numxvec(this.a * 1.0/wvl * DwvDxi, wv);
        let DCstrvDxj = numxvec(this.a * 1.0/wvl * DwvDxj, wv);
        let DCstrvDxk = numxvec(this.a * 1.0/wvl * DwvDxk, wv);

        let DCshDxi = sumvecs(numxvec(this.a * DwuDxi, wv), numxvec(this.a * DwvDxi, wu));
        let DCshDxj = sumvecs(numxvec(this.a * DwuDxj, wv), numxvec(this.a * DwvDxj, wu));
        let DCshDxk = sumvecs(numxvec(this.a * DwuDxk, wv), numxvec(this.a * DwvDxk, wu));

        let dotCstru = dot(DCstruDxi, this.V[i]) + dot(DCstruDxj, this.V[j]) + dot(DCstruDxk, this.V[k]);
        let dotCstrv = dot(DCstrvDxi, this.V[i]) + dot(DCstrvDxj, this.V[j]) + dot(DCstrvDxk, this.V[k]);
        let dotCsh = dot(DCshDxi, this.V[i]) + dot(DCshDxj, this.V[j]) + dot(DCshDxk, this.V[k]);

        let fi = sumarrayofvecs([numxvec(- ks * Cstru, DCstruDxi), numxvec(-ks * Cstrv, DCstrvDxi), numxvec( -ks * Csh, DCshDxi), 
            numxvec(- kd * dotCstru, DCstruDxi), numxvec(- kd * dotCstrv, DCstrvDxi), numxvec(-kd * dotCsh, DCshDxi)]); //damping contribution
        let fj = sumarrayofvecs([numxvec(- ks * Cstru, DCstruDxj), numxvec(-ks * Cstrv, DCstrvDxj), numxvec( -ks * Csh, DCshDxj), 
            numxvec(- kd * dotCstru, DCstruDxj), numxvec(- kd * dotCstrv, DCstrvDxj), numxvec(-kd * dotCsh, DCshDxj)]);
        let fk = sumarrayofvecs([numxvec(- ks * Cstru, DCstruDxk), numxvec(-ks * Cstrv, DCstrvDxk), numxvec( -ks * Csh, DCshDxk), 
            numxvec(- kd * dotCstru, DCstruDxk), numxvec(- kd * dotCstrv, DCstrvDxk), numxvec(-kd * dotCsh, DCshDxk)]);
        
        this.F[i] = sumvecs(this.F[i], fi);
        this.F[j] = sumvecs(this.F[j], fj);
        this.F[k] = sumvecs(this.F[k], fk);

    }

    /**
     * Считаем все силы.
    */

    calcForces()
    {
        for (let vertex = 0; vertex < this.F.length; vertex++)
            this.F[vertex] = sumvecs(this.F[vertex], numxvec(this.g, [0, 0, 1]));

        for (let triangle = 0; triangle < this.triangles.length; triangle++)
            this.calcTriangle(triangle);

        for (let pair = 0; pair < this.tri_pairs.length; pair++)
            this.calcTriPair(pair);
    }

    /**
     * Обновляем содержание вершинного объекта, который используется при рендеринге.
     */
    update()
    {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        gl.bufferData(gl.ARRAY_BUFFER, this.X, gl.STREAM_DRAW);
    }

    /**
     * Процедура симуляции.
     * @param {Float} dt 
     */

    simulate(dt)
    {
        let b = [];

        this.zero();

        this.calcForces();

        for (let k = 0; k < this.F.length; k++)
            b.push(numxvec((this.minv * dt), this.F[k]));
        
        let dv =  bdmxvec(this.S, b);

        let X0, X1, X2;

        for (let k = 0; k < this.V.length; k++)
        {
            /**
             * Интегрирование Верле.
             */
            X0 = this.X[5*k];
            X1 = this.X[5*k + 1];
            X2 = this.X[5*k + 2];
            
            this.X[5*k] = 2 * this.X[5*k] - this.Y[3*k] + dt *dt * dv[k][0];
            this.X[5*k + 1] = 2 * this.X[5*k + 1] - this.Y[3*k + 1] + dt * dt * dv[k][1];
            this.X[5*k + 2] = 2 * this.X[5*k + 2] - this.Y[3*k + 2] + dt * dt * dv[k][2];

            
            //Так как используем демпинг, должны сохранять данные о скорости.

            this.V[k] = numxvec(1/(1 * dt), difvecs([this.X[5*k], this.X[5*k + 1], this.X[5*k + 2]], [X0, X1, X2]));


            this.Y[3*k] = X0;
            this.Y[3*k + 1] = X1;
            this.Y[3*k + 2] = X2;

        }
    }



    /**
     * Процедура отрисовки.
     * @param {WebGLContext} gl 
     */
    drawCloth(gl)
    {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.ebo);
        gl.drawElements(gl.TRIANGLES, this.tri_inds.length, gl.UNSIGNED_SHORT, 0);
    }
}