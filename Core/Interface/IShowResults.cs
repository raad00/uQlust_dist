using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;

namespace uQlustCore.Interface
{
    public interface IShowResults
    {
        void Show(List<KeyValuePair<string,DataTable>> T);
        void ShowException(Exception ex);
    }
}
