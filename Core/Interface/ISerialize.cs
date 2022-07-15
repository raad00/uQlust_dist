using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace uQlustCore.Interface
{
    public interface ISerialize
    {
        ClusterOutput outCl { get; set;}
        void ISaveBinary(string fileName);
        ISerialize ILoadBinary(string fileName);
        
        List<string> ITest(string fileName);
    }
}
