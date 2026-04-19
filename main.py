from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from numbers import Number
from loguru import logger
import argparse


logger.add("filter_fastq.log", level="INFO", 
           format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
           rotation="10 MB", retention="1 month")


def filter_fastq(input_fastq: str, output_fastq: str = 'output.fastq', 
                 gc_bounds: Union[int, float, tuple] = (0, 100), 
                 length_bounds: Union[int, tuple] = (0, 2**32), 
                 quality_threshold: int = 0):
    """
    Filter sequences from a FASTQ file by GC content, length, and quality.
    
    Reads records from `input_fastq` and writes those meeting the provided
    criteria to `output_fastq` in FASTQ format.
    
    Args:
        input_fastq (str): Path to input FASTQ file.
        output_fastq (str): Path where filtered records will be written.
        gc_bounds (Union[int, float, tuple]): GC percentage boundaries. If a
            number, lower bound is 0; if a tuple, it should be (min, max).
        length_bounds (Union[int, tuple]): Sequence length boundaries. If an
            int, lower bound is 0; if a tuple, it should be (min, max).
        quality_threshold (int): Minimum average PHRED quality score.
    """

    def parse_bounds(bounds):
        if isinstance(bounds, Number):
            return 0, bounds
        return bounds[0], bounds[1]
    

    lower_gc_bound, upper_gc_bound = parse_bounds(gc_bounds)
    lower_length_bound, upper_length_bound = parse_bounds(length_bounds)


    total_records = 0
    passed_records = 0
          
    with open(input_fastq) as input_handle, \
         open(output_fastq, 'w') as output_handle:
        
        for seq_record in SeqIO.parse(input_handle, "fastq"):
            total_records += 1
            
            gc_content = gc_fraction(seq_record.seq) * 100
            seq_length = len(seq_record)
            qualities = seq_record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities)
            
            if (lower_gc_bound <= gc_content <= upper_gc_bound and
                lower_length_bound <= seq_length <= upper_length_bound and
                quality_threshold <= avg_quality):
                SeqIO.write(seq_record, output_handle, "fastq")
                passed_records += 1
    
    logger.info(f"Total processed records: {total_records}. \n \
                Matching the filtering criteria: {passed_records}.")


def parse_bounds_from_cli(bounds_values, default_min, default_max, param_name):
    """
    Parse command line bounds arguments.
    
    Args:
        bounds_values: List of values from command line (or None)
        default_min: Default minimum value
        default_max: Default maximum value
        param_name: Name of parameter for error messages
    
    Returns:
        Tuple of (min, max)
    """
    # Если аргумент не указан, возвращаем значения по умолчанию
    if bounds_values is None or len(bounds_values) == 0:
        return default_min, default_max
    
    # Если указано одно число
    if len(bounds_values) == 1:
        return default_min, bounds_values[0]
    
    # Если указано два числа
    if len(bounds_values) == 2:
        min_val, max_val = bounds_values
        if min_val > max_val:
            raise ValueError(f"Error in {param_name}: min value ({min_val}) cannot be greater than max value ({max_val})")
        return min_val, max_val
    
    # Если указано больше двух чисел
    raise ValueError(f"Error in {param_name}: expected 0, 1 or 2 values, got {len(bounds_values)}")


def create_parser():

    parser = argparse.ArgumentParser(
                        description='Filter FASTQ sequences by GC content, length, and quality')

    parser.add_argument('input_fastq', type=str, help='Path to input FASTQ file')
    parser.add_argument('-o','--output', type=str, default='output.fastq', dest='output_fastq', 
                        help='Path to the output FASTQ with filtered fastq-sequences. \
                        Default: `output.fastq`.')

    filter_group = parser.add_argument_group('FASTQ-filtration')
    filter_group.add_argument('-g', '--gc', nargs="*", type=int, help='GC percentage boundaries. \
                        If one number, lower bound is 0; if two numbers, it must be `min` `max`. \
                        Default: `0` `100`.')
    filter_group.add_argument('-l', '--length', nargs="*", type=int, help='Sequence length boundaries. \
                        If one number, lower bound is 0; if two numbers, it must be `min` `max`. \
                        Default: `0` `2**32`.')
    filter_group.add_argument('-q', '--quality', type=int, default=0, help='Minimum average PHRED quality score. \
                        Default: 0.')

    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    gc_min, gc_max = parse_bounds_from_cli(
        args.gc, 
        default_min=0, 
        default_max=100, 
        param_name="GC bounds"
    )

    length_min, length_max = parse_bounds_from_cli(
        args.length, 
        default_min=0, 
        default_max=2**32, 
        param_name="Length bounds"
    )


    params = params = {
        'gc_lower_bound': gc_min,
        'gc_upper_bound': gc_max,
        'length_lower_bound': length_min,
        'length_upper_bound': length_max,
        'quality': args.quality
    }

    try:
        negative_params = [p for p, value in params.items() if value < 0]
        if negative_params:
            raise ValueError(f"Negative values found in: {', '.join(negative_params)}.")
    except ValueError as error:
        logger.error(f"{str(error)} The filtering parameters cannot be negative")
        raise


    filter_fastq(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=(gc_min, gc_max),
        length_bounds=(length_min, length_max),
        quality_threshold=args.quality
    )
    
