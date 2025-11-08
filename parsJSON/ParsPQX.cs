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
                PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
                PropertyNameCaseInsensitive = true
            };
            return JsonSerializer.Deserialize<PQX_model>(jsonString, options);;
        }
    }
}
