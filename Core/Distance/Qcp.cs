using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace uQlustCore.Distance
{
    class Qcp
    {
        //public  double[,] rot = new double[3, 3];
        public double rot00, rot01, rot02, rot10, rot11, rot12, rot20, rot21, rot22; //instead of array, because of linux!
        //double[] C = new double[4];
        double C0, C1, C2;
        //float[,] covMat;
        int size;
        double E0;
        double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
        double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
        double mxEigenV;
        public Qcp(float[,] covMat, double E0, int size)
        {
            /*this.covMat = new double[covMat.GetLength(0), covMat.GetLength(1)];
            for (int i = 0; i < covMat.GetLength(0); i++)
                for (int j = 0; j < covMat.GetLength(1); j++)
                    this.covMat[i, j] = (double)covMat[i, j];*/

            Sxx = covMat[0, 0]; Sxy = covMat[0, 1]; Sxz = covMat[0, 2];
            Syx = covMat[1, 0]; Syy = covMat[1, 1]; Syz = covMat[1, 2];
            Szx = covMat[2, 0]; Szy = covMat[2, 1]; Szz = covMat[2, 2];

            this.size = size;
            this.E0 = E0;
        }



        void PrepCalc()
        {
            /*Sxx = covMat[0, 0]; Sxy = covMat[0, 1]; Sxz = covMat[0, 2];
            Syx = covMat[1, 0]; Syy = covMat[1, 1]; Syz = covMat[1, 2];
            Szx = covMat[2, 0]; Szy = covMat[2, 1]; Szz = covMat[2, 2];*/

            Sxx2 = Sxx * Sxx;
            Syy2 = Syy * Syy;
            Szz2 = Szz * Szz;

            Sxy2 = Sxy * Sxy;
            Syz2 = Syz * Syz;
            Sxz2 = Sxz * Sxz;

            Syx2 = Syx * Syx;
            Szy2 = Szy * Szy;
            Szx2 = Szx * Szx;

            SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
            Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

            C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
            C1 = 8.0 * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz);

            SxzpSzx = Sxz + Szx;
            SyzpSzy = Syz + Szy;
            SxypSyx = Sxy + Syx;
            SyzmSzy = Syz - Szy;
            SxzmSzx = Sxz - Szx;
            SxymSyx = Sxy - Syx;
            SxxpSyy = Sxx + Syy;
            SxxmSyy = Sxx - Syy;
            Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

            C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
                 + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
                + (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz)) * (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
                + (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz)) * (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz)) * (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz)) * (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));


        }

        public double CalcRmsd()
        {
            int i;

            double oldg = 0.0;
            double b, a, x2, delta;
            double evalprec = 1e-11;

            PrepCalc();

            mxEigenV = E0;
            for (i = 0; i < 50; ++i)
            {
                oldg = mxEigenV;
                x2 = mxEigenV * mxEigenV;
                b = (x2 + C2) * mxEigenV;
                a = b + C1;
                delta = ((a * mxEigenV + C0) / (2.0 * x2 * mxEigenV + b + a));
                mxEigenV -= delta;
                if (Math.Abs(mxEigenV - oldg) < Math.Abs(evalprec * mxEigenV))
                    break;
            }

            if (i == 20)
            {
//Console.WriteLine("\nMore than " + i + " iterations needed!");
   //             return -1;
            }

            /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
            return Math.Sqrt(Math.Abs(2.0 * (E0 - mxEigenV) / size));
        }
        public bool CalcRotMatrix()
        {
            double qsqr;
            double q1, q2, q3, q4, normq;
            double a11, a12, a13, a14, a21, a22, a23, a24;
            double a31, a32, a33, a34, a41, a42, a43, a44;
            double a2, x2, y2, z2;
            double xy, az, zx, ay, yz, ax;
            double a3344_4334, a3244_4234, a3243_4233, a3143_4133, a3144_4134, a3142_4132;
            double evecprec = 1e-6;

            //PrepCalc();

            rot00 = rot11 = rot22 = 1.0;
            rot01 = rot02 = rot10 = rot12 = rot20 = rot21 = 0.0;

            double rmsd=CalcRmsd();
            if(rmsd<0)
                return false;
               
            a11 = SxxpSyy + Szz - mxEigenV; a12 = SyzmSzy; a13 = -SxzmSzx; a14 = SxymSyx;
            a21 = SyzmSzy; a22 = SxxmSyy - Szz - mxEigenV; a23 = SxypSyx; a24 = SxzpSzx;
            a31 = a13; a32 = a23; a33 = Syy - Sxx - Szz - mxEigenV; a34 = SyzpSzy;
            a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
            a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44 - a42 * a34;
            a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43 - a41 * a33;
            a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42 - a41 * a32;
            q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
            q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
            q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
            q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

            qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

            /* The following code tries to calculate another column in the adjoint matrix when the norm of the 
               current column is too small.
               Usually this commented block will never be activated.  To be absolutely safe this should be
               uncommented, but it is most likely unnecessary.  
            */
            if (qsqr < evecprec)
            {
                q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
                q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
                q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
                q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
                qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                if (qsqr < evecprec)
                {
                    double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                    double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                    double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                    q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                    q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                    q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                    q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                    if (qsqr < evecprec)
                    {
                        q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                        q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                        q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                        q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                        qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                        if (qsqr < evecprec)
                        {
                            /* if qsqr is still too small, return the identity matrix. */
                            return false;
                        }
                    }
                }
            }

            normq = Math.Sqrt(qsqr);
            q1 /= normq;
            q2 /= normq;
            q3 /= normq;
            q4 /= normq;

            a2 = q1 * q1;
            x2 = q2 * q2;
            y2 = q3 * q3;
            z2 = q4 * q4;

            xy = q2 * q3;
            az = q1 * q4;
            zx = q4 * q2;
            ay = q1 * q3;
            yz = q3 * q4;
            ax = q1 * q2;


            rot00 = a2 + x2 - y2 - z2;
            rot01 = 2 * (xy + az);
            rot02 = 2 * (zx - ay);
            rot10 = 2 * (xy - az);
            rot11 = a2 - x2 + y2 - z2;
            rot12 = 2 * (yz + ax);
            rot20 = 2 * (zx + ay);
            rot21 = 2 * (yz - ax);
            rot22 = a2 - x2 - y2 + z2;
            return true;

        }

        
    }
}
