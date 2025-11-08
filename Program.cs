using celiang;
using LL.parsJSON;
namespace LL
{
    class Chengxu
    {
        public static void Main(String[] args)
        {
            //var pqx = ParsPQX.ParseJson(@"D:\73\function\hdm\bin\Release\net8.0-windows\publish\2d\lnglat\广州.pqx");
            var rd = celiang.LL.pqx;

           

            // 查找名为"江灵南"的道路
            Road jiangLingNanRoad = rd.roads?.FirstOrDefault(road => road.name == "江灵南");

            if (jiangLingNanRoad?.pqx != null)
            {
                // 此时 jiangLingNanRoad.pqx 已经是 List<List<double>> 类型
                // 如果需要转换为二维数组 double[][]
                double[][] twoDimensionalArray = jiangLingNanRoad.pqx
                    .Select(innerList => innerList.ToArray())
                    .ToArray();

                // 或者直接使用 List<List<double>>
                List<List<double>> pqxList = jiangLingNanRoad.pqx;
            }
            Console.ReadKey();
        }
    }
}