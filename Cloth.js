class Cloth
{
    constructor(size, gl)
    {
        // a = 1.0 / (2 * (size * size));
        // g = 0.98;
        // ks = 100;
        // kd = 10;
        // bu = 0.9;
        // bv = 0.9;
        // kb = 10;
        //a = 1.0 / (2 * (size * size));

        // a = 1.0;
        // g = 9.8;
        // ks = 1.4;
        // kd = 0.05;
        // kb = 3;
        // kdb = 0.5;
        // bu = 0.9;
        // bv = 0.9;
        // //n = 7
        this.size = size;
        this.gl = gl;
        this.X = new Float32Array(size * size * 5);
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

        this.minv = 1
        this.a = 1.0;
        this.g = 9.8;
        this.ks = 1000.0;
        this.kd = 15.8;
        this.kb = 50.0;
        this.kdb = 5.8;
        this.bu = 1;
        this.bv = 1;
        
        for (let i = 0; i < size; i++)
        {
            for (let j = 0; j < size; j++)
            {
                this.u.push([20 * j/(size - 1), 20 * i/(size - 1)]);
                this.X[5*(i * size + j)] = 20 * j/(size - 1);
                this.X[5*(i * size + j) + 1] = 20 * i/(size - 1);
                this.X[5*(i * size + j) + 2] = 0;
                this.X[5*(i * size + j) + 3] = 20 * j/(size - 1);
                this.X[5*(i * size + j) + 4] = 20 * i/(size - 1);

                this.F.push([0, 0, 0]);
                this.V.push([0, 0, 0]);
                this.z.push([0, 0, 0]);
                
                if (i == 0)
                    this.S.push([[0, 0, 0], [0, 0, 0], [0, 0, 0]]); //diag(1.0) - glm::vec3(0, 0, 1) * glm::vec3(0, 0, 1));
                else
                    this.S.push([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
            }                

        }
        
            
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

        let NTR = 2 * (size - 1); //number of triangles in a row
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

        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo); 
        gl.bufferData(gl.ARRAY_BUFFER, this.X, gl.STREAM_DRAW);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.ebo);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.tri_inds), gl.STATIC_DRAW);

        gl.enableVertexAttribArray(0);
        gl.enableVertexAttribArray(1);
        gl.vertexAttribPointer(0, 3, gl.FLOAT, false, this.X.BYTES_PER_ELEMENT * 5, 0);
        gl.vertexAttribPointer(1, 2, gl.FLOAT, false, this.X.BYTES_PER_ELEMENT * 5, this.X.BYTES_PER_ELEMENT * 3);

    }

    zero()
    {
        let n;
        
        n = this.F.length;
        for (let k = 0; k < n; k++)
            this.F[k] = [0, 0, 0];

    }

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


    calcTriPair(pair)
    {
        let i, j, k, m, l;
        let inters = [];
        let ex = [];

        let kb = this.kb;
        let kdb = this.kdb;

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

        let dx1 = [this.X[j*5] - this.X[i*5], this.X[j*5 + 1] - this.X[i*5 + 1], this.X[j*5 + 2] - this.X[i*5 + 2]];
        let dx2 = [this.X[k*5] - this.X[i*5], this.X[k*5 + 1] - this.X[i*5 + 1], this.X[k*5 + 2] - this.X[i*5 + 2]];
        let dx3 = [this.X[k*5] - this.X[j*5], this.X[k*5 + 1] - this.X[j*5 + 1], this.X[k*5 + 2] - this.X[j*5 + 2]];
        let dx4 = [this.X[m*5] - this.X[j*5], this.X[m*5 + 1] - this.X[j*5 + 1], this.X[m*5 + 2] - this.X[j*5 + 2]];

        let N1 = cross(dx1, dx2);
        let N2 = cross(dx3, dx4);
        if (dot(N1, N2) < 0)
            N2 = numxvec(-1, N2);
        
        let Cb = 1 - (1.0/length(N1)) * (1.0/length(N2)) * dot(N1, N2);           

        let Dn1Dxi = [0, difvecs(dx1, dx2)[2], difvecs(dx2, dx1)[1], difvecs(dx2, dx1)[2], 0, difvecs(dx1, dx2)[0], difvecs(dx1, dx2)[1], difvecs(dx2, dx1)[0], 0];
        let Dn1Dxj = [0, dx2[2], -dx2[1], -dx2[2], 0, dx2[0], dx2[1], -dx2[0], 0];
        let Dn1Dxk = [0, -dx1[2], dx1[1], dx1[2], 0, -dx1[0], -dx1[1], dx1[0], 0];

        let Dn2Dxj = [0, difvecs(dx3, dx4)[2], difvecs(dx4, dx3)[1], difvecs(dx4, dx3)[2], 0, difvecs(dx3, dx4)[0], difvecs(dx3, dx4)[1], difvecs(dx4, dx3)[0], 0];
        let Dn2Dxk = [0, dx4[2], -dx4[1], -dx4[2], 0, dx4[0], dx4[1], -dx4[0], 0];
        let Dn2Dxm = [0, -dx3[2], dx3[1], dx3[2], 0, -dx3[0], -dx3[1], dx3[0], 0];
        
        let DCbDxi = numxvec((-1.0/length(N1)) * (1.0/length(N2)), matxvec(Dn1Dxi, N2));
        let DCbDxj = numxvec((-1.0/length(N1)) * (1.0/length(N2)), sumvecs(matxvec(Dn1Dxj, N2), matxvec(Dn2Dxj, N1)));
        let DCbDxk = numxvec((-1.0/length(N1)) * (1.0/length(N2)), sumvecs(matxvec(Dn1Dxk, N2), matxvec(Dn2Dxk, N1)));
        let DCbDxm = numxvec((-1.0/length(N1)) * (1.0/length(N2)), matxvec(Dn2Dxm, N1));


        let dotCb = dot(DCbDxi, this.V[i]) + dot(DCbDxj, this.V[j]) + dot(DCbDxk, this.V[k]) + dot(DCbDxm, this.V[m]);

        //calculate contribution to forces from triangles
        let fi = sumvecs(numxvec(-kb * Cb, DCbDxi), numxvec(-kdb * dotCb, DCbDxi)); //damping contribution
        let fj = sumvecs(numxvec(-kb * Cb, DCbDxj), numxvec(-kdb * dotCb, DCbDxj));
        let fk = sumvecs(numxvec(-kb * Cb, DCbDxk), numxvec(-kdb * dotCb, DCbDxk));
        let fm = sumvecs(numxvec(-kb * Cb, DCbDxm), numxvec(-kdb * dotCb, DCbDxm));
        
        this.F[i] = sumvecs(this.F[i], fi);
        this.F[j] = sumvecs(this.F[j], fj);
        this.F[k] = sumvecs(this.F[k], fk);
        this.F[m] = sumvecs(this.F[m], fm);

    }


    calcTriangle(triangle)
    {
        let i = this.triangles[triangle][0];
        let j = this.triangles[triangle][1];
        let k = this.triangles[triangle][2];

        let ks = this.ks;
        let kd = this.kd;

        let dx1 = [this.X[5*j] - this.X[5*i], this.X[5*j + 1] - this.X[5*i + 1], this.X[5*j + 2] - this.X[5*i + 2]];
        let dx2 = [this.X[5*k] - this.X[5*i], this.X[5*k + 1] - this.X[5*i + 1], this.X[5*k + 2] - this.X[5*i + 2]];

        let wu = sumvecs(numxvec(this.Uinv[triangle][0], dx1), numxvec(this.Uinv[triangle][2], dx2));
        let wv = sumvecs(numxvec(this.Uinv[triangle][1], dx1), numxvec(this.Uinv[triangle][3], dx2));
        let wul = length(wu);
        let wvl = length(wv);

        let Cstru = this.a * (wul - this.bu);
        let Cstrv = this.a * (wvl - this.bv);
        let Csh =  this.a * dot(wu, wv); 

        let DwuDxi = diag(- this.Uinv[triangle][0] - this.Uinv[triangle][2]);
        let DwuDxj = diag(this.Uinv[triangle][0]);
        let DwuDxk = diag(this.Uinv[triangle][2]);
        let DwvDxi = diag(-this.Uinv[triangle][1] - this.Uinv[triangle][3]);
        let DwvDxj = diag(this.Uinv[triangle][1]);
        let DwvDxk = diag(this.Uinv[triangle][3]);

        let DCstruDxi = numxvec(this.a * 1.0/wul, matxvec(DwuDxi, wu));
        let DCstruDxj = numxvec(this.a * 1.0/wul, matxvec(DwuDxj, wu));
        let DCstruDxk = numxvec(this.a * 1.0/wul, matxvec(DwuDxk, wu));

        let DCstrvDxi = numxvec(this.a * 1.0/wvl, matxvec(DwvDxi,wv));
        let DCstrvDxj = numxvec(this.a * 1.0/wvl, matxvec(DwvDxj,wv));
        let DCstrvDxk = numxvec(this.a * 1.0/wvl, matxvec(DwvDxk,wv));

        let DCshDxi = sumvecs(numxvec(this.a, matxvec(DwuDxi, wv)), numxvec(this.a, matxvec(DwvDxi, wu)));
        let DCshDxj = sumvecs(numxvec(this.a, matxvec(DwuDxj, wv)), numxvec(this.a, matxvec(DwvDxj, wu)));
        let DCshDxk = sumvecs(numxvec(this.a, matxvec(DwuDxk, wv)), numxvec(this.a, matxvec(DwvDxk, wu)));

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

    calcForces()
    {
        for (let vertex = 0; vertex < this.F.length; vertex++)
            this.F[vertex] = sumvecs(this.F[vertex], numxvec(this.g, [0, 0, 1]));

        for (let triangle = 0; triangle < this.triangles.length; triangle++)
            this.calcTriangle(triangle);

        for (let pair = 0; pair < this.tri_pairs.length; pair++)
            this.calcTriPair(pair);
    }

    update()
    {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        gl.bufferData(gl.ARRAY_BUFFER, this.X, gl.STREAM_DRAW);
    }

    simulate(dt)
    {
        let b = [];

        this.zero();

        this.calcForces();

        for (let k = 0; k < this.F.length; k++)
            b.push(numxvec((this.minv * dt), this.F[k]));
        
        let dv =  bdmxvec(this.S, b);

        for (let k = 0; k < this.V.length; k++)
            this.V[k] = sumvecs(this.V[k],dv[k]);

        for (let k = 0; k < this.V.length; k++)
        {
            this.X[5*k] = this.X[5*k] + dt * (this.V[k][0]);
            this.X[5*k + 1] = this.X[5*k + 1] + dt * (this.V[k][1]);
            this.X[5*k + 2] = this.X[5*k + 2] + dt * (this.V[k][2]);
        }
    }




    drawCloth(gl)
    {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.ebo);
        gl.drawElements(gl.TRIANGLES, this.tri_inds.length, gl.UNSIGNED_SHORT, 0);
    }
}