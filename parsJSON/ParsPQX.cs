using System.Text.Json;
namespace LL.parsJSON
{
    public class ParsPQX
    {
        public static PQX_model? ParseJson(string filePath)
        {
            string jsonString = File.ReadAllText(filePath);
            var options = new JsonSerializerOptions
            {
                PropertyNamingPolicy = JsonNamingPolicy.CamelCase, // 属性名驼峰式匹配
                PropertyNameCaseInsensitive = true, // 忽略大小写
                DefaultIgnoreCondition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull // 忽略空值
            };
            return JsonSerializer.Deserialize<PQX_model>(jsonString, options);;
        }
    }
}
