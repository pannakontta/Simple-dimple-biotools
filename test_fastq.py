import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Импортируем тестируемые функции
from main import parse_bounds_from_cli, filter_fastq


class TestParseBoundsFromCLI:
    
    def test_no_bounds_returns_defaults(self):
        result = parse_bounds_from_cli(None, 0, 100, "test")
        assert result == (0, 100)
        
        result = parse_bounds_from_cli([], 0, 100, "test")
        assert result == (0, 100)
    
    def test_single_bound_returns_default_min(self):
        result = parse_bounds_from_cli([50], 0, 100, "test")
        assert result == (0, 50)
    
    def test_two_bounds_valid(self):
        result = parse_bounds_from_cli([20, 80], 0, 100, "test")
        assert result == (20, 80)
    
    def test_min_bound_greater_than_max_raises_error(self):
        with pytest.raises(ValueError, match=r"Error in test bounds: min value \(\d+\) cannot be greater than max value \(\d+\)"):
            parse_bounds_from_cli([100, 50], 0, 100, "test bounds")
    
    def test_too_many_bounds_raises_error(self):
        with pytest.raises(ValueError, match="expected 0, 1 or 2 values, got 3"):
            parse_bounds_from_cli([1, 2, 3], 0, 100, "test")


class TestFilterFastqFileOperations:
    
    @pytest.fixture
    def sample_fastq_content(self):
        records = [
            SeqRecord(Seq("AAAAA"), id="seq1", description="", 
                     letter_annotations={"phred_quality": [30, 30, 30, 30, 30]}),
            SeqRecord(Seq("GCGCG"), id="seq2", description="",
                     letter_annotations={"phred_quality": [40, 40, 40, 40, 40]}),
            SeqRecord(Seq("ATATAT"), id="seq3", description="",
                     letter_annotations={"phred_quality": [10, 10, 10, 10, 10, 10]}),
        ]
        return records
    
    def test_filter_fastq_input_and_output(self, sample_fastq_content, tmp_path):
        """Тест: чтение из входного файла и запись в выходной (чтение/запись файла)"""
        
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"
        
        with open(input_file, 'w') as f:
            SeqIO.write(sample_fastq_content, f, "fastq")
        
        filter_fastq(str(input_file), str(output_file), 
                    gc_bounds=(0, 100), length_bounds=(0, 100), quality_threshold=0)
        
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
        with open(output_file) as f:
            filtered_records = list(SeqIO.parse(f, "fastq"))
            assert len(filtered_records) == 3
    
    def test_filter_fastq_with_criteria_filters_correctly(self, sample_fastq_content, tmp_path):
        
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"
        
        with open(input_file, 'w') as f:
            SeqIO.write(sample_fastq_content, f, "fastq")
        
        filter_fastq(str(input_file), str(output_file), 
                    gc_bounds=(50, 100), length_bounds=(0, 100), quality_threshold=0)
        
        with open(output_file) as f:
            filtered_records = list(SeqIO.parse(f, "fastq"))
            assert len(filtered_records) == 1
            assert filtered_records[0].id == "seq2"
    
    def test_filter_fastq_empty_output_when_no_matches(self, sample_fastq_content, tmp_path):
        
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"
        
        with open(input_file, 'w') as f:
            SeqIO.write(sample_fastq_content, f, "fastq")
        
        # Нереальные требования
        filter_fastq(str(input_file), str(output_file), 
                    gc_bounds=(1000, 2000), length_bounds=(0, 100), quality_threshold=0)
        
        with open(output_file) as f:
            content = f.read()
            assert content == "" or len(list(SeqIO.parse(content, "fastq"))) == 0

