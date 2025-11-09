using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;
using celiang;

namespace LL
{
    internal class JiaoDian
    {


        public double[,] JD2PQX(double[,] data)
        {
            List<double> listXY1 = [];
            var lr18 = data;
            var sk = lr18[0, 2];
            var JD2Mileage = 0.0;
            double hzx = 0, hzy = 0;
            (double X, double Y) JD1;
            (double X, double Y) JD2;
            (double X, double Y) JD3;
            for (int i = 1; i < lr18.GetLength(0) - 1; i++)
            {
                JD1 = (lr18[i - 1, 0], lr18[i - 1, 1]);

                if (i > 1)
                {
                    JD1 = (hzx, hzy);
                }
                JD2 = (lr18[i, 0], lr18[i, 1]);

                JD3 = (lr18[i + 1, 0], lr18[i + 1, 1]);

                var R = lr18[i, 4];
                var Ls1 = lr18[i, 2];
                var Ls2 = lr18[i, 3];
                var xyzy = 1;
                var azimuth12 = calculateAzimuth(JD1, JD2);//方位角
                var azimuth23 = calculateAzimuth(JD2, JD3);
                var alpha = azimuth23 - azimuth12;
                if (alpha < 0)
                {
                    alpha = -alpha;
                }
                if (alpha > Math.PI)
                {
                    alpha = Math.PI * 2 - alpha;
                }
                var area = calculateTriangleArea(JD1.X, JD1.Y, JD2.X, JD2.Y, JD3.X, JD3.Y);
                if (area < 0)
                    xyzy = -1;
                var (T1, T2, Ly, L, a2, a3, a4, a5, a6) = calculateCurveElements(Ls1, Ls2, R, alpha);//切线
                double ZH, HY, QZ, YH, HZ;
                double zhx, zhy;
                var dist1 = calculateDistance(JD1.X, JD1.Y, JD2.X, JD2.Y);
                JD2Mileage = sk + dist1;
                ZH = JD2Mileage - T1;
                HY = ZH + Ls1;
                QZ = ZH + Ls1 + (L - Ls2 - Ls1) / 2;
                YH = ZH + L - Ls2;
                HZ = ZH + L;

                //"jd" + i + "前直线"
                if (dist1 > T1 && dist1 - T1 > 0.01)
                {
                    double xycd = dist1 - T1;

                    listXY1.AddRange([sk, JD1.X, JD1.Y, azimuth12, xycd, 0, 0, 0]);
                }
                zhx = JD2.X - T1 * Math.Cos(azimuth12);
                zhy = JD2.Y - T1 * Math.Sin(azimuth12);
                var zhk = sk + dist1 - T1;
                // "jd" + i + "一缓", 
                if (Ls1 != 0)
                {
                    listXY1.AddRange([zhk, zhx, zhy, azimuth12, Ls1, 0, R, xyzy]);
                }
                //var xy = {
                //    x: 0,
                //    y: 0,
                //    "fwj": 0,
                //    rad: 0
                //}
            ;
                //"jd" + i + "圆弧", 
                double[] xy = [0, 0, 0];
                if (Ly != 0 && Ls1 != 0)
                {
                    xy = celiang.LL.Zs(ZH, zhx, zhy, azimuth12, Ls1, 0, R, xyzy, ZH + Ls1, 0, 90);
                    //[targetX, targetY, targetAzimuth];
                    listXY1.AddRange([HY, xy[0], xy[1], xy[2], Ly, R, R, xyzy]);
                }
                else
                {
                    //"圆弧jd" + i,
                    //xy = {
                    //x: zhx,
                    //        y: zhy,
                    //        "fwj": 0,
                    //        rad: azimuth12
                    //        }
                    ;
                    listXY1.AddRange([HY, zhx, zhy, azimuth12, Ly, R, R, xyzy]);
                }
                //"jd" + i + "二缓", 
                if (Ls2 != 0)
                {
                    xy = celiang.LL.Zs(HY, xy[0], xy[1], xy[2], Ly, R, R, xyzy, HY + Ly, 0, 90);
                    listXY1.AddRange([YH, xy[0], xy[1], xy[2], Ls2, R, 0, xyzy]);
                }
                hzx = JD2.X + T2 * Math.Cos(azimuth23);
                hzy = JD2.Y + T2 * Math.Sin(azimuth23);
                sk = zhk + L;
                //"jd" + i + "终点的直线", 
                if (i == lr18.GetLength(0) - 2)
                {
                    var dist2 = calculateDistance(JD2.X, JD2.Y, JD3.X, JD3.Y);
                    if (dist2 - T2 > 0.01)
                    {
                        listXY1.AddRange([sk, hzx, hzy, azimuth23, dist2 - T2, 0, 0, 0]);
                    }
                }
            }
            int rows = listXY1.Count / 8;
            int cols = 8;
            double[,] pqx = new double[rows,8];
            int lie = 0;
           for (int i = 0; i < listXY1.Count; i ++)
            {
                int row = i / cols;  // 计算行索引： 
                int col = i % cols;  // 计算列索引： 
                pqx[row, col] = listXY1[i];
            }
            for (int i = 0; i < pqx.GetLength(0); i++)
            {
                pqx[i, 3] = celiang.LL.radiansToDMS(pqx[i, 3]);
            }
            return pqx;
        }



        public double calculateAzimuth((double X, double Y) point1, (double X, double Y) point2)
        {
            var dx = point2.X - point1.X;
            var dy = point2.Y - point1.Y;
            var azimuth = Math.Atan2(dy, dx);
            if (azimuth < 0)
            {
                azimuth += 2 * Math.PI;
            }
            return azimuth;
        }

        public double calculateTriangleArea(double x1, double y1, double x2, double y2, double x3, double y3)
        {

            return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
        }

        public (double, double, double, double, double, double, double, double, double) calculateCurveElements(double Ls1, double Ls2, double R, double alpha)
        {
            double T1 = 0;
            double T2 = 0;
            double E, L, J;
            var m1 = Ls1 / 2 - Math.Pow(Ls1, 3) / (240 * Math.Pow(R, 2)) - Math.Pow(Ls1, 5) / (34560 * Math.Pow(R, 4));
            var m2 = Ls2 / 2 - Math.Pow(Ls2, 3) / (240 * Math.Pow(R, 2)) - Math.Pow(Ls2, 5) / (34560 * Math.Pow(R, 4));
            var p1 = Math.Pow(Ls1, 2) / (24 * R) - Math.Pow(Ls1, 4) / (2688 * R * R * R);
            var p2 = Math.Pow(Ls2, 2) / (24 * R) - Math.Pow(Ls2, 4) / (2688 * R * R * R);
            T1 = m1 + (R + p2 - (R + p1) * Math.Cos(alpha)) / Math.Sin(alpha);
            T2 = m2 + (R + p1 - (R + p2) * Math.Cos(alpha)) / Math.Sin(alpha);
            E = (R + (p1 + p2) / 2) / Math.Cos(alpha / 2) - R;
            var beta01 = Ls1 / (2 * R);
            var beta02 = Ls2 / (2 * R);
            double Ly = R * (alpha - beta01 - beta02);
            L = Ls1 + Ls2 + Ly;
            J = T1 + T2 - L;
            return (T1, T2, Ly, L, Ly, p1, p2, beta01, beta02);
        }
        public double calculateDistance(double x1, double y1, double x2, double y2)
        {
            var dx = Math.Pow(x2 - x1, 2);
            var dy = Math.Pow(y2 - y1, 2);
            return Math.Sqrt(dx + dy);
        }

    }
}