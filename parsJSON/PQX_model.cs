namespace LL.parsJSON
{
    public class PQX_model
    {
        public string proj { get; set; }
        public int center { get; set; }
        public bool south { get; set; }
        public List<DataItem> data { get; set; }
        public List<Road> roads { get; set; }
    }

    public class DataItem
    {
        public double B { get; set; }
        public double L { get; set; }
        public double x { get; set; }
        public double y { get; set; }
    }

    public class Road
    {
        public string name { get; set; }
        public List<List<double>> DL { get; set; }
        public List<List<double>> jd { get; set; }
        public List<List<double>> sqx { get; set; }
        public List<List<double>> pqx { get; set; }
    }

}
