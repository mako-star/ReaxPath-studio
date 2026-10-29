"""ReaxFF bond trajectory analysis inspired by ChemTraYzer's workflow."""

from reaxnet.analysis import AnalysisResult, AnalyzerConfig, analyze_file, analyze_frames
from reaxnet.streaming import StreamingAnalysisResult, analyze_file_streaming

__all__ = [
    "AnalysisResult",
    "AnalyzerConfig",
    "StreamingAnalysisResult",
    "analyze_file",
    "analyze_file_streaming",
    "analyze_frames",
]
__version__ = "0.1.0"
