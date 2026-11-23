package picard.flow;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

public class FlowBasedReadTpDecodeTest {

    // ticket example: QUAL DDDD:::*:D:D::, mq:D
    private static final String SEQ = "CAGAGCCACCCACG";
    private static final String QUAL = "DDDD:::*:D:D::";

    private static SAMRecord makeRecord(final String mq, final String qual, final Object tp) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        if (mq != null) {
            rg.setAttribute(FlowBasedRead.MAX_QUAL_READ_GROUP_TAG, mq);
        }
        header.addReadGroup(rg);
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName("r1");
        rec.setReadUnmappedFlag(true);
        rec.setReadString(SEQ.substring(0, qual.length()));
        rec.setBaseQualityString(qual);
        rec.setAttribute("RG", "rg1");
        if (tp != null) {
            rec.setAttribute("tp", tp);
        }
        return rec;
    }

    @Test
    public void testOldFormat() {
        final byte[] dense = {0, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, 0, 1, 1};
        Assert.assertEquals(FlowBasedRead.decodeTp(makeRecord(null, QUAL, dense)), dense);
    }

    @Test
    public void testSparse() {
        Assert.assertEquals(FlowBasedRead.decodeTp(makeRecord("D", QUAL, "+--+++++")),
                new byte[]{0, 0, 0, 0, 1, -1, -1, 1, 1, 0, 1, 0, 1, 1});
    }

    @Test
    public void testAbsentTagWithMq() {
        Assert.assertEquals(FlowBasedRead.decodeTp(makeRecord("D", "DDDD", null)), new byte[4]);
    }

    @Test
    public void testEmptyString() {
        Assert.assertEquals(FlowBasedRead.decodeTp(makeRecord("D", "DDDD", "")), new byte[4]);
    }

    @DataProvider
    public Object[][] invalid() {
        return new Object[][]{
                {"D", QUAL, "+--+++"},       // too short
                {"D", QUAL, "+--++++++"},    // leftover
                {"D", QUAL, "+--++0++"},     // bad char
                {"D", "DDDD", "+"},          // leftover on all-max read
                {null, QUAL, "+--+++++"},    // no mq
                {null, QUAL, null},          // no mq, no tp
                {"D", QUAL, new byte[3]},    // old format, length mismatch
        };
    }

    @Test(dataProvider = "invalid", expectedExceptions = PicardException.class)
    public void testInvalid(final String mq, final String qual, final Object tp) {
        FlowBasedRead.decodeTp(makeRecord(mq, qual, tp));
    }
}
