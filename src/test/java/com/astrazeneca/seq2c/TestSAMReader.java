package com.astrazeneca.seq2c;


import htsjdk.samtools.*;
import org.mockito.invocation.InvocationOnMock;
import org.mockito.stubbing.Answer;
import org.powermock.api.mockito.PowerMockito;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class TestSAMReader implements SamReader {

    @Override
    public SAMFileHeader getFileHeader() {
        final SAMSequenceRecord mockedRecord = mock(SAMSequenceRecord.class);
        when(mockedRecord.getSequenceName()).thenReturn("chr1");

        SAMSequenceDictionary mockedDictionary = mock(SAMSequenceDictionary.class);
        when(mockedDictionary.getSequences()).thenReturn(new ArrayList<SAMSequenceRecord>() {{
            add(mockedRecord);
        }});

        SAMFileHeader mockHeader = mock(SAMFileHeader.class);
        when(mockHeader.getSequenceDictionary()).thenReturn(mockedDictionary);

        return mockHeader;
    }

    @Override
    public Type type() {
        return null;
    }

    @Override
    public String getResourceDescription() {
        return null;
    }

    @Override
    public boolean hasIndex() {
        return false;
    }

    @Override
    public Indexing indexing() {
        return null;
    }

    @Override
    public SAMRecordIterator iterator() {
        return null;
    }

    @Override
    public SAMRecordIterator query(String sequence, int start, int end, boolean contained) {
        final List<SAMRecord> records = new ArrayList<SAMRecord>(){{

            add(createRecord(20,1, new ArrayList<CigarElement>(){{
                add(new CigarElement(2, CigarOperator.M));
                add(new CigarElement(3, CigarOperator.D));
            }}));
            add(createRecord(15, 1, new ArrayList<CigarElement>(){{
                add(new CigarElement(2, CigarOperator.S));
                add(new CigarElement(2, CigarOperator.M));
                add(new CigarElement(1, CigarOperator.D));
            }}));
            add(createRecord(15, 255, new ArrayList<CigarElement>(){{
                add(new CigarElement(2, CigarOperator.M));
                add(new CigarElement(1, CigarOperator.D));
                add(new CigarElement(2, CigarOperator.S));
            }}));
        }};

        final Iterator<SAMRecord> recordIterator = records.iterator();
        SAMRecordIterator iterator = PowerMockito.mock(SAMRecordIterator.class);
        when(iterator.hasNext()).thenAnswer(new Answer<Object>() {
            @Override
            public Object answer(InvocationOnMock invocation) throws Throwable {
                return recordIterator.hasNext();
            }
        });
        when(iterator.next()).thenAnswer(new Answer<Object>() {
            @Override
            public Object answer(InvocationOnMock invocation) throws Throwable {
                return recordIterator.next();
            }
        });

        return iterator;
    }

    static SAMRecord createRecord(int alignmentStart, int flag, List<CigarElement> cigarElements) {
        SAMRecord record = mock(SAMRecord.class);
        when(record.getFlags()).thenReturn(flag);
        when(record.getAlignmentStart()).thenReturn(alignmentStart);
        when(record.getReferenceName()).thenReturn("test_name");
        when(record.getMateReferenceName()).thenReturn("test_name");
        when(record.getCigar()).thenReturn(new Cigar(cigarElements));
        return record;
    }

    @Override
    public SAMRecordIterator queryOverlapping(String sequence, int start, int end) {
        return null;
    }

    @Override
    public SAMRecordIterator queryContained(String sequence, int start, int end) {
        return null;
    }

    @Override
    public SAMRecordIterator query(QueryInterval[] intervals, boolean contained) {
        return null;
    }

    @Override
    public SAMRecordIterator queryOverlapping(QueryInterval[] intervals) {
        return null;
    }

    @Override
    public SAMRecordIterator queryContained(QueryInterval[] intervals) {
        return null;
    }

    @Override
    public SAMRecordIterator queryUnmapped() {
        return null;
    }

    @Override
    public SAMRecordIterator queryAlignmentStart(String sequence, int start) {
        return null;
    }

    @Override
    public SAMRecord queryMate(SAMRecord rec) {
        return null;
    }

    @Override
    public void close() throws IOException {
    }
}
