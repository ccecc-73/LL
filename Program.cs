using celiang;
using LL.parsJSON;
namespace LL
{
    class Chengxu
    {
        public static void Main(String[] args)
        {
            ////var pqx = ParsPQX.ParseJson(@"D:\73\function\hdm\bin\Release\net8.0-windows\publish\2d\lnglat\广州.pqx");
            //var rd = celiang.LL.pqx;



            //// 查找名为"江灵南"的道路
            //Road jiangLingNanRoad = rd.roads?.FirstOrDefault(road => road.name == "江灵南");

            //if (jiangLingNanRoad?.pqx != null)
            //{
            //    // 此时 jiangLingNanRoad.pqx 已经是 List<List<double>> 类型
            //    // 如果需要转换为二维数组 double[][]
            //    double[][] twoDimensionalArray = jiangLingNanRoad.pqx
            //        .Select(innerList => innerList.ToArray())
            //        .ToArray();

            //    // 或者直接使用 List<List<double>>
            //    List<List<double>> pqxList = jiangLingNanRoad.pqx;
            //}


            double[] source = [1710.9090, 884.1963, 1717.1560, 892.0050, 1708.0775, 898.0678];
            double[] target = [1728.5216, 878.5348, 1737.1560, 892.0050, 1721.6961, 900.1340];

            double[] cs = celiang.LL.Cs4(source, target);
            celiang.LL.FourParameterTransform(1710.9090, 884.1963, cs[0], cs[1], cs[2], cs[3]);


            for (int i = 0;i<target.Length;i+=2)
            {
                var res = celiang.LL.FourParameterTransform(source[i], source[i + 1], cs[0], cs[1], cs[2], cs[3]);
                Console.WriteLine($"源点({source[i]:F3}, {source[i + 1]:F3}) => 目标点({res[0]:F3}, {res[1]:F3})， 期望目标点({target[i]:F3}, {target[i + 1]:F3})");
            }
            Console.ReadKey();
        }
    }
}