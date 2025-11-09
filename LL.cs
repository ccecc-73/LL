using LL.parsJSON;

namespace celiang
{
    /// <summary>
    /// 
    /// </summary>
    public class LL
    {
        public static PQX_model? pqx { get; set; }

        static LL()
        {
           // pqx = ParsPQX.ParseJson(@"D:\73\function\hdm\bin\Release\net8.0-windows\publish\2d\lnglat\广州.pqx");
        }
        //1
        #region   数字转桩号格式 
        /// <summary>
        /// 将数字转换为 K 格式，例如：12345 -> K12+345
        /// </summary>
        /// <param name="meters">输入的数字</param>
        /// <returns>格式化后的字符串</returns>
        /// <example>LL.Num2K(12345)</example>
        public string Num2K(double meters)
        {
            int km = (int)Math.Floor(meters / 1000);
            double m = meters - km * 1000;
            m = Math.Round(m, 2);

            if (m == (int)m)
            {
                // 整数米，如 2.00 → k1+002
                return $"K{km}+{(int)m:D3}";
            }
            else
            {
                // 小数米，如 2.35 → k0+002.35
                return $"K{km}+{m:000.00}";
            }
        }
        #endregion

        //2
        #region 格式为 DD.MMSS转弧度
        /// <summary>
        /// 将度分秒（DMS）格式转换为弧度
        /// 假设输入的 dms 是一个浮点数，格式为 DD.MMSS
        /// </summary>
        /// <param name="dms">度分秒格式的数值</param>
        /// <returns>对应的弧度值</returns>
        internal static double DmsToRadians(double dms)
        {
            int degrees = (int)dms;
            double fractional = dms - degrees;
            int minutes = (int)(fractional * 100);
            double seconds = (fractional * 100 - minutes) * 100;

            double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
            return totalDegrees * Math.PI / 180;
        }
        #endregion

        //3
        #region  距离和方位角
        /// <summary>
        /// 计算两点之间的极坐标距离和角度（弧度）
        /// </summary>
        /// <param name="x0">起点 X 坐标</param>
        /// <param name="y0">起点 Y 坐标</param>
        /// <param name="x1">终点 X 坐标</param>
        /// <param name="y1">终点 Y 坐标</param>
        /// <returns>包含两个元素的数组：[距离(cd), 角度(hd]) ]</returns>
        internal static double[] Fwj(double x0, double y0, double x1, double y1)
        {
            double x = x1 - x0;
            double y = y1 - y0;
            double cd = Math.Sqrt(x * x + y * y); // 距离
            double hd = Math.Atan2(y, x);         // 角度（弧度）

            if (hd < 0)
            {
                hd += 2 * Math.PI; // 确保角度在 [0, 2π) 范围内
            }
            return [cd, hd];
        }
        #endregion

        //4
        #region 由桩号计算xy
        /// <summary>
        /// k2xy
        /// </summary>
        /// <param name="xyk">起始弧长</param>
        /// <param name="xyx">起始X坐标</param>
        /// <param name="xyy">起始Y坐标</param>
        /// <param name="xyhd">起始方位角（弧度）</param>
        /// <param name="xycd">当前圆弧半径</param>
        /// <param name="xyqdr">起始圆弧半径</param>
        /// <param name="xyzdr">结束圆弧半径</param>
        /// <param name="xyzy">方向系数</param>
        /// <param name="jsk">目标弧长</param>
        /// <param name="jsb">偏移距离</param>
        /// <param name="jd">偏移角度（度）</param>
        /// <returns>包含目标X, 目标Y, 目标方位角（弧度）的对象，同时附加属性 x, y, fwj（度）, rad（弧度）</returns>
        internal static double[] Zs(double xyk, double xyx, double xyy, double xyhd, double xycd, double xyqdr, double xyzdr, double xyzy, double jsk, double jsb, double jd)
        {
            // 第一个条件分支  圆弧
            if (Math.Abs(xyqdr - xyzdr) < 0.01 && xyqdr != 0)
            {
                double centerX = xyx + xyqdr * Math.Cos(xyhd + xyzy * Math.PI / 2);
                double centerY = xyy + xyqdr * Math.Sin(xyhd + xyzy * Math.PI / 2);
                double deltaArcLength = jsk - xyk;
                double deltaAngle = deltaArcLength / xyqdr * xyzy;
                double targetAzimuth = xyhd + deltaAngle;
                if (targetAzimuth < 0)
                    targetAzimuth += 2 * Math.PI;
                double angle = xyhd + xyzy * Math.PI / 2 + deltaAngle;
                double targetX = centerX - xyqdr * Math.Cos(angle) + jsb * Math.Cos(angle - xyzy * Math.PI / 2 + jd * Math.PI / 180);
                double targetY = centerY - xyqdr * Math.Sin(angle) + jsb * Math.Sin(angle - xyzy * Math.PI / 2 + jd * Math.PI / 180);

                return [targetX, targetY, targetAzimuth];
            }

            // 第二个条件分支 直线
            if (xyqdr < 0.01 && xyzdr < 0.01 && xyzy < 0.01)
            {
                double targetX = xyx + (jsk - xyk) * Math.Cos(xyhd) + jsb * Math.Cos(xyhd + jd * Math.PI / 180);
                double targetY = xyy + (jsk - xyk) * Math.Sin(xyhd) + jsb * Math.Sin(xyhd + jd * Math.PI / 180);

                return [targetX, targetY, xyhd];
            }

            // 第三个条件分支（默认情况） 缓和曲线
            if (xyqdr < 0.001) xyqdr = 99999999;
            if (xyzdr < 0.001) xyzdr = 99999999;

            double f0 = xyhd;
            double q = xyzy;
            double c = 1 / xyqdr;
            double d = (xyqdr - xyzdr) / (2 * xycd * xyqdr * xyzdr);

            // 初始化 rr 和 vv 数组，索引从1开始，因此大小为5（索引1-4）
            double[] rr = new double[5];
            double[] vv = new double[5];

            rr[1] = 0.1739274226;
            rr[2] = 0.3260725774;
            rr[3] = rr[2];
            rr[4] = rr[1];

            vv[1] = 0.0694318442;
            vv[2] = 0.3300094782;
            vv[3] = 1 - vv[2];
            vv[4] = 1 - vv[1];

            double w = jsk - xyk;
            double xs = 0;
            double ys = 0;

            for (int i = 1; i < 5; i++)
            {
                double ff = f0 + q * vv[i] * w * (c + vv[i] * w * d);
                xs += rr[i] * Math.Cos(ff);
                ys += rr[i] * Math.Sin(ff);
            }

            double fhz3 = f0 + q * w * (c + w * d);
            if (fhz3 < 0)
                fhz3 += 2 * Math.PI;
            if (fhz3 >= 2 * Math.PI)
                fhz3 -= 2 * Math.PI;

            double fhz1 = xyx + w * xs + jsb * Math.Cos(fhz3 + jd * Math.PI / 180);
            double fhz2 = xyy + w * ys + jsb * Math.Sin(fhz3 + jd * Math.PI / 180);
            return [fhz1, fhz2, fhz3];
        }
        #endregion

        //5
        #region 单条路线单个坐标计算
        /// <summary>
        /// 计算指定里程处的坐标和方位角
        /// </summary>
        /// <param name="pqx">路径点数组，每个元素是一个包含多个属性的数组</param>
        /// <param name="k">目标里程</param>
        /// <param name="b">宽度参数</param>
        /// <param name="z">角度参数</param>
        /// <returns>包含目标X坐标、Y坐标和方位角的数组</returns>
        public static double[] Dantiaoxianludange(double[,] pqx, double k, double b, double z)
        {
            for (int i = 0; i < pqx.GetLength(0); i++)
            {
                double dtk = pqx[i, 0];
                double dtx = pqx[i, 1];
                double dty = pqx[i, 2];
                double dtfwj = pqx[i, 3];
                double dtcd = pqx[i, 4];
                double dtr1 = pqx[i, 5];
                double dtr2 = pqx[i, 6];
                double dtzy = pqx[i, 7];
                if (k >= dtk && k <= dtk + dtcd)
                {
                    double hudu = DmsToRadians(dtfwj);
                    double[] jsxy1 = Zs(dtk, dtx, dty, hudu, dtcd, dtr1, dtr2, dtzy, k, b, z);
                    return [Math.Round(jsxy1[0], 3), Math.Round(jsxy1[1], 3), jsxy1[2]];
                }
            }
            return [0, 0, 0];
        }
        #endregion

        //6
        #region xy=>k
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pqx"></param>
        /// <param name="fsx"></param>
        /// <param name="fsy"></param>
        /// <returns></returns>
        public static double[] Fs(double fsx, double fsy, double[,] pqx)
        {
            dynamic jljd = Fwj(pqx[0, 1], pqx[0, 2], fsx, fsy);
            double k = pqx[0, 0];
            double hudu = DmsToRadians(pqx[0, 3]);
            double cz = jljd[0] * Math.Cos(jljd[1] - hudu);
            double pj = jljd[0] * Math.Sin(jljd[1] - hudu);
            int hang = pqx.GetLength(0) - 1;
            double qdlc = pqx[0, 0];
            double zdlc = pqx[hang, 0] + pqx[hang, 4];
            int jisuancishu = 0;
            while (Math.Abs(cz) > 0.01)
            {
                k += cz;
                jisuancishu += 1;
                if (k < qdlc)
                {
                    return [-1, -1];
                }
                if (k > zdlc)
                {
                    return [-2, -2];
                }
                if (jisuancishu > 15)
                {
                    return [-3, -3];
                }
                double[] xy = Dantiaoxianludange(pqx, k, 0, 0);
                jljd = Fwj((double)xy[0], (double)xy[1], fsx, fsy);
                cz = jljd[0] * Math.Cos(jljd[1] - xy[2]);
                pj = jljd[0] * Math.Sin(jljd[1] - xy[2]);
            }
            return [Math.Round(k, 3), Math.Round(pj, 3)];
        }
        #endregion

        #region 
        private static double Gaocheng(double bpdlc, double bpdgc, double r, double qp, double hp, double t, double k)
        {
            double f = qp - hp;
            r = r * Math.Abs(f) / f;
            double x;
            if (k <= bpdlc - t)
            {
                x = 0;
            }
            else if (k >= bpdlc + t)
            {
                x = 0;
                qp = hp;
            }
            else
            {
                x = k - bpdlc + t;
            }

            return (bpdgc - (bpdlc - k) * qp - Math.Pow(x, 2) / 2 / r);
        }
        #endregion

        //7高程计算
        #region h
        /// <summary>
        /// 
        /// </summary>
        /// <param name="k"></param>
        /// <param name="sqxb"></param>
        /// <returns></returns>
        public static double H(double k, double[,] sqxb)
        {
            double hp = 0;
            int length = sqxb.GetLength(0);
            for (int i = 1; i < length - 1; i++)
            {
                double r = sqxb[i, 2];
                if (r < 0.001)
                    r = 0.001;
                double qp = (sqxb[i, 1] - sqxb[i - 1, 1]) / (sqxb[i, 0] - sqxb[i - 1, 0]);
                hp = (sqxb[i + 1, 1] - sqxb[i, 1]) / (sqxb[i + 1, 0] - sqxb[i, 0]);
                double f = qp - hp;
                double t = r * Math.Abs(f) / 2;
                if (k <= sqxb[i, 0] + t)
                    return Math.Round(Gaocheng(sqxb[i, 0], sqxb[i, 1], r, qp, hp, t, k), 3);
            }
            //the last
            if (k <= sqxb[length - 1, 0])
            {
                return Math.Round(sqxb[length - 1, 1] + (k - sqxb[length - 1, 0]) * hp, 3);
            }
            return -1;
        }
        #endregion
        //高斯投影
        /// <summary>
        /// 
        /// </summary>
        /// <param name="L">经度</param>
        /// <param name="B">纬度</param>
        /// <param name="lonCenter">中心经度,为空的时候自动计算</param>
        /// <returns>[北坐标,东坐标,中心经度]</returns>
        public static double[] Gauss_proj(double L, double B, double lonCenter = 360)
        {
            double pi = 3.141592653589793238463;
            double p0 = 206264.8062470963551564;
            double e = 0.00669438002290;
            double e1 = 0.00673949677548;
            double b = 6356752.3141;
            double a = 6378137.0;
            B = B * pi / 180;
            L = L * pi / 180;
            double L_num;
            double L_center;
            if (lonCenter == 360)
            {
                L_num = Math.Floor(L * 180 / pi / 3.0 + 0.5);
                L_center = 3 * L_num;
            }
            else
            {
                L_center = (double)lonCenter;
            }
            double l = (L / pi * 180 - L_center) * 3600;
            double M0 = a * (1 - e);
            double M2 = 3.0 / 2.0 * e * M0;
            double M4 = 5.0 / 4.0 * e * M2;
            double M6 = 7.0 / 6.0 * e * M4;
            double M8 = 9.0 / 8.0 * e * M6;
            double a0 = M0 + M2 / 2.0 + 3.0 / 8.0 * M4 + 5.0 / 16.0 * M6 + 35.0 / 128.0 * M8;
            double a2 = M2 / 2.0 + M4 / 2 + 15.0 / 32.0 * M6 + 7.0 / 16.0 * M8;
            double a4 = M4 / 8.0 + 3.0 / 16.0 * M6 + 7.0 / 32.0 * M8;
            double a6 = M6 / 32.0 + M8 / 16.0;
            double a8 = M8 / 128.0;
            double Xz = a0 * B - a2 / 2.0 * Math.Sin(2 * B) + a4 / 4.0 * Math.Sin(4 * B) - a6 / 6.0 * Math.Sin(6 * B) + a8 / 8.0 * Math.Sin(8 * B);
            double c = a * a / b;
            double V = Math.Sqrt(1 + e1 * Math.Cos(B) * Math.Cos(B));
            double N = c / V;
            double t = Math.Tan(B);
            double n = e1 * Math.Cos(B) * Math.Cos(B);
            double m1 = N * Math.Cos(B);
            double m2 = N / 2.0 * Math.Sin(B) * Math.Cos(B);
            double m3 = N / 6.0 * Math.Pow(Math.Cos(B), 3) * (1 - t * t + n);
            double m4 = N / 24.0 * Math.Sin(B) * Math.Pow(Math.Cos(B), 3) * (5 - t * t + 9 * n);
            double m5 = N / 120.0 * Math.Pow(Math.Cos(B), 5) * (5 - 18 * t * t + Math.Pow(t, 4) + 14 * n - 58 * n * t * t);
            double m6 = N / 720.0 * Math.Sin(B) * Math.Pow(Math.Cos(B), 5) * (61 - 58 * t * t + Math.Pow(t, 4));
            double x = Xz + m2 * l * l / Math.Pow(p0, 2) + m4 * Math.Pow(l, 4) / Math.Pow(p0, 4) + m6 * Math.Pow(l, 6) / Math.Pow(p0, 6);
            double y0 = m1 * l / p0 + m3 * Math.Pow(l, 3) / Math.Pow(p0, 3) + m5 * Math.Pow(l, 5) / Math.Pow(p0, 5);
            double y = y0 + 500000;
            return [x, y, L_center];
        }
        //高斯 反投影
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x">北坐标</param>
        /// <param name="y">东坐标</param>
        /// <param name="l0">中心经度</param>
        /// <returns>[经度,纬度]</returns>
        public static double[] Gauss_unproj(double x, double y, double l0)
        {
            double pi = 3.141592653589793238463;
            double e = 0.00669438002290;
            double e1 = 0.00673949677548;
            double b = 6356752.3141;
            double a = 6378137.0;
            double y1 = y - 500000;
            double M0 = a * (1 - e);
            double M2 = 3.0 / 2.0 * e * M0;
            double M4 = 5.0 / 4.0 * e * M2;
            double M6 = 7.0 / 6.0 * e * M4;
            double M8 = 9.0 / 8.0 * e * M6;
            double a0 = M0 + M2 / 2.0 + 3.0 / 8.0 * M4 + 5.0 / 16.0 * M6 + 35.0 / 128.0 * M8;
            double a2 = M2 / 2.0 + M4 / 2 + 15.0 / 32.0 * M6 + 7.0 / 16.0 * M8;
            double a4 = M4 / 8.0 + 3.0 / 16.0 * M6 + 7.0 / 32.0 * M8;
            double a6 = M6 / 32.0 + M8 / 16.0;
            double Bf = x / a0;
            double B0 = Bf;
            while (Math.Abs(Bf - B0) > 0.0000001 || B0 == Bf)
            {
                B0 = Bf;
                double FBf = -a2 / 2.0 * Math.Sin(2 * B0) + a4 / 4.0 * Math.Sin(4 * B0) - a6 / 6.0 * Math.Sin(6 * B0);
                Bf = (x - FBf) / a0;
            }
            double t = Math.Tan(Bf);
            double c = a * a / b;
            double V = Math.Sqrt(1 + e1 * Math.Cos(Bf) * Math.Cos(Bf));
            double N = c / V;
            double M = c / Math.Pow(V, 3);
            double n = e1 * Math.Cos(Bf) * Math.Cos(Bf);
            double n1 = 1 / (N * Math.Cos(Bf));
            double n2 = -t / (2.0 * M * N);
            double n3 = -(1 + 2 * t * t + n) / (6.0 * Math.Pow(N, 3) * Math.Cos(Bf));
            double n4 = t * (5 + 3 * t * t + n - 9 * n * t * t) / (24.0 * M * Math.Pow(N, 3));
            double n5 = (5 + 28 * t * t + 24 * Math.Pow(t, 4) + 6 * n + 8 * n * t * t) / (120.0 * Math.Pow(N, 5) * Math.Cos(Bf));
            double n6 = -t * (61 + 90 * t * t + 45 * Math.Pow(t, 4)) / (720.0 * M * Math.Pow(N, 5));
            double B = (Bf + n2 * y1 * y1 + n4 * Math.Pow(y1, 4) + n6 * Math.Pow(y1, 6)) / pi * 180;
            double L0 = l0;
            double l = n1 * y1 + n3 * Math.Pow(y1, 3) + n5 * Math.Pow(y1, 5);
            double L = L0 + l / pi * 180;
            return [L, B];
        }

        // 经纬度转UTM坐标
        /// <summary>
        /// 
        /// </summary>
        /// <param name="longitude"></param>
        /// <param name="latitude"></param>
        /// <returns></returns>
        public static double[] Utm_proj(double longitude, double latitude)
        {
            double EQUATORIAL_RADIUS = 6378137.0;
            double FLATTENING = 1 / 298.257223563;
            double ECC_SQUARED = 2 * FLATTENING - Math.Pow(FLATTENING, 2);
            double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
            double SCALE_FACTOR = 0.9996;
            double FALSE_EASTING = 500000.0;
            double FALSE_NORTHING_S = 10000000.0;
            int zoneNumber = (int)Math.Floor((longitude + 180) / 6) + 1;
            double centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;
            double latRad = (latitude) * Math.PI / 180.0;
            double lonRad = (longitude) * Math.PI / 180.0;
            double lonCenterRad = (centralMeridian) * Math.PI / 180.0;
            double N = EQUATORIAL_RADIUS / Math.Sqrt(1 - ECC_SQUARED * Math.Pow(Math.Sin(latRad), 2));
            double T = Math.Pow(Math.Tan(latRad), 2);
            double C = ECC_PRIME_SQUARED * Math.Pow(Math.Cos(latRad), 2);
            double A = (lonRad - lonCenterRad) * Math.Cos(latRad);
            double M = EQUATORIAL_RADIUS * ((1 - ECC_SQUARED / 4 - 3 * Math.Pow(ECC_SQUARED, 2) / 64 - 5 * Math.Pow(ECC_SQUARED, 3) / 256) * latRad - (3 * ECC_SQUARED / 8 + 3 * Math.Pow(ECC_SQUARED, 2) / 32 + 45 * Math.Pow(ECC_SQUARED, 3) / 1024) * Math.Sin(2 * latRad) + (15 * Math.Pow(ECC_SQUARED, 2) / 256 + 45 * Math.Pow(ECC_SQUARED, 3) / 1024) * Math.Sin(4 * latRad) - (35 * Math.Pow(ECC_SQUARED, 3) / 3072) * Math.Sin(6 * latRad));
            double easting = SCALE_FACTOR * N * (A + (1 - T + C) * Math.Pow(A, 3) / 6 + (5 - 18 * T + Math.Pow(T, 2) + 72 * C - 58 * ECC_PRIME_SQUARED) * Math.Pow(A, 5) / 120) + FALSE_EASTING;
            double northing = SCALE_FACTOR * (M + N * Math.Tan(latRad) * (Math.Pow(A, 2) / 2 + (5 - T + 9 * C + 4 * Math.Pow(C, 2)) * Math.Pow(A, 4) / 24 + (61 - 58 * T + Math.Pow(T, 2) + 600 * C - 330 * ECC_PRIME_SQUARED) * Math.Pow(A, 6) / 720));
            if (latitude < 0) northing += FALSE_NORTHING_S;
            return [northing, easting, zoneNumber];
        }

        // UTM坐标转经纬度
        /// <summary>
        /// 
        /// </summary>
        /// <param name="northing"></param>
        /// <param name="easting"></param>
        /// <param name="zoneNumber"></param>
        /// <param name="isNorthern"></param>
        /// <returns></returns>
        public static double[] Utm_unproj(double northing, double easting, int zoneNumber, bool isNorthern)
        {
            double EQUATORIAL_RADIUS = 6378137.0;
            double FLATTENING = 1 / 298.257223563;
            double ECC_SQUARED = 2 * FLATTENING - Math.Pow(FLATTENING, 2);
            double ECC_PRIME_SQUARED = ECC_SQUARED / (1 - ECC_SQUARED);
            double SCALE_FACTOR = 0.9996;
            double FALSE_EASTING = 500000.0;
            double FALSE_NORTHING_S = 10000000.0;
            double x = easting - FALSE_EASTING;
            double y = isNorthern ? northing : northing - FALSE_NORTHING_S;
            double centralMeridian = (zoneNumber - 1) * 6 - 180 + 3;
            double lonCenterRad = (centralMeridian) * Math.PI / 180.0;
            double M = y / SCALE_FACTOR;
            double mu = M / (EQUATORIAL_RADIUS * (1 - ECC_SQUARED / 4 - 3 * Math.Pow(ECC_SQUARED, 2) / 64
                - 5 * Math.Pow(ECC_SQUARED, 3) / 256));
            double e1 = (1 - Math.Sqrt(1 - ECC_SQUARED)) / (1 + Math.Sqrt(1 - ECC_SQUARED));
            double phi1Rad = mu + (3 * e1 / 2 - 27 * Math.Pow(e1, 3) / 32) * Math.Sin(2 * mu)
                + (21 * Math.Pow(e1, 2) / 16 - 55 * Math.Pow(e1, 4) / 32) * Math.Sin(4 * mu)
                + (151 * Math.Pow(e1, 3) / 96) * Math.Sin(6 * mu);
            double N1 = EQUATORIAL_RADIUS / Math.Sqrt(1 - ECC_SQUARED * Math.Pow(Math.Sin(phi1Rad), 2));
            double T1 = Math.Pow(Math.Tan(phi1Rad), 2);
            double C1 = ECC_PRIME_SQUARED * Math.Pow(Math.Cos(phi1Rad), 2);
            double R1 = EQUATORIAL_RADIUS * (1 - ECC_SQUARED) / Math.Pow(1 - ECC_SQUARED * Math.Pow(Math.Sin(phi1Rad), 2), 1.5);
            double D = x / (N1 * SCALE_FACTOR);
            double latRad = phi1Rad - (N1 * Math.Tan(phi1Rad) / R1)
                * (Math.Pow(D, 2) / 2 - (5 + 3 * T1 + 10 * C1 - 4 * Math.Pow(C1, 2) - 9 * ECC_PRIME_SQUARED)
                * Math.Pow(D, 4) / 24 + (61 + 90 * T1 + 298 * C1 + 45 * Math.Pow(T1, 2)
                - 252 * ECC_PRIME_SQUARED - 3 * Math.Pow(C1, 2)) * Math.Pow(D, 6) / 720);
            double lonRad = lonCenterRad + (D - (1 + 2 * T1 + C1) * Math.Pow(D, 3) / 6
                + (5 - 2 * C1 + 28 * T1 - 3 * Math.Pow(C1, 2) + 8 * ECC_PRIME_SQUARED + 24 * Math.Pow(T1, 2))
                * Math.Pow(D, 5) / 120) / Math.Cos(phi1Rad);
            return [(latRad) * 180 / Math.PI, (lonRad) * 180 / Math.PI];
        }

        //四参数
        // 应用转换公式
        //double convertedX = params.scale* (x* Math.cos(params.rotation) - y* Math.sin(params.rotation)) + params.deltaX;
        //        double convertedY = params.scale* (x* Math.sin(params.rotation) + y* Math.cos(params.rotation)) + params.deltaY;
        /// <summary>
        /// 计算两组坐标系下的四参数转换参数
        /// </summary>
        /// <param name="source">[x1,y1,x2,y2...]</param>
        /// <param name="target">[x1,y1,x2,y2...]</param>
        /// <returns>[deltaX, deltaY, rotation(弧度), scale]</returns>
        /// <exception cref="ArithmeticException"></exception>
        public static double[] Cs4(double[] source, double[] target)
        {
            if (source == null || target == null || source.Length != target.Length)
            {
                Console.WriteLine("坐标数组长度必须相等");
                return [0, 0, 0, 1];
            }
            if (source.Length < 4 || source.Length % 2 != 0)
            {
                Console.WriteLine("至少需要2个点且坐标为偶数");
                return [0, 0, 0, 1];
            }
            int pointCount = source.Length / 2;
            double sumX1 = 0, sumY1 = 0, sumX2 = 0, sumY2 = 0;
            for (int i = 0; i < pointCount; i++)
            {
                sumX1 += source[2 * i];
                sumY1 += source[2 * i + 1];
                sumX2 += target[2 * i];
                sumY2 += target[2 * i + 1];
            }
            double meanX1 = sumX1 / pointCount;
            double meanY1 = sumY1 / pointCount;
            double meanX2 = sumX2 / pointCount;
            double meanY2 = sumY2 / pointCount;
            double[] centeredSource = new double[source.Length];
            double[] centeredTarget = new double[target.Length];
            for (int i = 0; i < pointCount; i++)
            {
                centeredSource[2 * i] = source[2 * i] - meanX1;
                centeredSource[2 * i + 1] = source[2 * i + 1] - meanY1;
                centeredTarget[2 * i] = target[2 * i] - meanX2;
                centeredTarget[2 * i + 1] = target[2 * i + 1] - meanY2;
            }
            double H11 = 0, H12 = 0, H21 = 0, H22 = 0;
            double B1 = 0, B2 = 0;
            for (int i = 0; i < pointCount; i++)
            {
                double x1 = centeredSource[2 * i];
                double y1 = centeredSource[2 * i + 1];
                double x2 = centeredTarget[2 * i];
                double y2 = centeredTarget[2 * i + 1];
                H11 += x1 * x1 + y1 * y1;
                H12 += 0;
                H21 += 0;
                H22 += x1 * x1 + y1 * y1;
                B1 += x1 * x2 + y1 * y2;
                B2 += x1 * y2 - y1 * x2;
            }
            double det = H11 * H22 - H12 * H21;
            if (Math.Abs(det) < 1e-15)
            {
                Console.WriteLine("矩阵奇异，无法求解参数");
            }
            double a = (H22 * B1 - H12 * B2) / det;
            double b = (-H21 * B1 + H11 * B2) / det;
            double scale = Math.Sqrt(a * a + b * b);
            double rotation = Math.Atan2(b, a);
            double deltaX = meanX2 - (a * meanX1 - b * meanY1);
            double deltaY = meanY2 - (b * meanX1 + a * meanY1);
             Console.WriteLine($"四参数计算:\ndeltaX: {deltaX}\ndeltaY:{deltaY}\nrotation (radians): {rotation}\nscale: {scale}");
            return [deltaX, deltaY, rotation, scale];
        }
        /// <summary>
        /// 已知四参数进行坐标转换
        /// </summary>
        /// <param name="x">原坐标x</param>
        /// <param name="y">原坐标y</param>
        /// <param name="deltaX">平移x</param>
        /// <param name="deltaY">平移y</param>
        /// <param name="rotation">旋转弧度</param>
        /// <param name="scale">缩放</param>
        /// <returns>[x,y]</returns>
        public static double[] FourParameterTransform(double x, double y, double deltaX, double deltaY, double rotation, double scale)
        {
            double convertedX = scale * (x * Math.Cos(rotation) - y * Math.Sin(rotation)) + deltaX;
            double convertedY = scale * (x * Math.Sin(rotation) + y * Math.Cos(rotation)) + deltaY;
            //Console.WriteLine($"转换后坐标:\nx: {convertedX}\ny: {convertedY}");
            return [convertedX, convertedY];
        }
    }
}