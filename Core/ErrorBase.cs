using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;

namespace uQlustCore
{
    public static class ErrorBase
    {
        private static List<string> errors = new List<string>();
        public static void ClearErrors()
        {
            errors.Clear();
        }
        public static void AddErrors(string error)
        {
            if(errors.Count<5000)
                errors.Add(error);
            if (errors.Count == 5000)
                errors.Add("There are much more errors but there is not enough room to store them");
        }
        public static List<string> GetErrors()
        {
            return errors;
        }
    }
}
