using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;


namespace uQlustCore
{
    /*public class GraphNode
    {
        public int x, y;
        public int areaLeft, areaRight;
        int square = 8;

        public bool MouseClick(int mouseX, int mouseY)
        {
            if (mouseX >= x - square && mouseX <= x  + square / 2 - square && mouseY >= y - square && mouseY <= y  + square / 2 - square)
                return true;
            return false;
        }

        public bool IsDrawAble(int recX, int recY, int width, int height)
        {
            if (x >= recX && x <= recX + width && y >= recY && y <= recY + height)
                return true;
            return false;
        }
        public void DrawNode(Graphics g,float lineThick)
        {
            Pen p;
            p = new Pen(Color.Black);
            p.Width = lineThick;              
            g.DrawRectangle(p, x - square, y - square, square, square);            
            g.DrawLine(p, x  + square / 2 - square, y  - square, x  + square / 2 - square, y + square - square);
            g.DrawLine(p, x  - square, y + square / 2 - square, x  + square - square, y + square / 2 - square);
        }
    }*/
}
