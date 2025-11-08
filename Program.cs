using celiang;
using LL.parsJSON;
namespace LL
{
    class Chengxu
    {
        public static void Main(String[] args)
        {
            var pqx = ParsPQX.ParseJson(@"D:\73\function\hdm\bin\Release\net8.0-windows\publish\2d\lnglat\广州.pqx");
            Console.WriteLine(value: pqx.proj);
            
            Console.ReadKey();
        }
    }
}