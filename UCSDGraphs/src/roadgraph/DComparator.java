package roadgraph;

import java.util.Comparator;

public  class DComparator implements Comparator<Double>
{
	

	@Override
	public int compare(Double x, Double y) {
		// TODO Auto-generated method stub
		if (x < y)
        {
            return -1;
        }
        if (x > y)
        {
            return 1;
        }
        return 0;
	}

	
}