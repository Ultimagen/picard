package picard.flow;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;

public class FlowBasedArgumentCollection implements Serializable {
    private static final long serialVersionUID = 0;

    public static final String FLOW_IGNORE_T0_TAG_LONG_NAME = "flow-ignore-t0-tag";
    public static final String FILLING_VALUE_LONG_NAME = "flow-fill-empty-bins-value";
    public static final String FLOW_KEEP_BOUNDARY_FLOWS_LONG_NAME = "keep-boundary-flows";
    private static final double DEFAULT_FILLING_VALUE = 0;
    private static final boolean DEFAULT_FLOW_IGNORE_T0_TAG = false;
    private static final boolean DEFAULT_FLOW_KEEP_BOUNDARY_FLOWS = false;

    @Advanced
    @Hidden
    @Argument(fullName = FLOW_IGNORE_T0_TAG_LONG_NAME, doc = "Ignore t0 tag in the read when create flow matrix (arcane/obsolete)", optional = true)
    public boolean ignoreT0Tag = DEFAULT_FLOW_IGNORE_T0_TAG;

    @Advanced
    @Hidden
    @Argument(fullName = FILLING_VALUE_LONG_NAME, doc = "Value to fill the zeros of the matrix with", optional=true)
    public double fillingValue = DEFAULT_FILLING_VALUE;

    @Advanced
    @Argument(fullName=FLOW_KEEP_BOUNDARY_FLOWS_LONG_NAME, doc="prevent spreading of boundary flows.", optional = true)
    public boolean keepBoundaryFlows = DEFAULT_FLOW_KEEP_BOUNDARY_FLOWS;

    public FlowBasedArgumentCollection() {}
}
