import os
import json
from pathlib import Path
import pdfplumber
import PyPDF2
import re
from urllib.parse import urlparse


def extract_pdf_text(pdf_path):
    """
    Extract text from PDF file preserving original formatting.
    Tries multiple methods to handle different PDF formats.

    Args:
        pdf_path (Path): Path to the PDF file

    Returns:
        str: Extracted text content
    """
    # Try pdfplumber first
    try:
        with pdfplumber.open(pdf_path) as pdf:
            text_content = ""
            for page in pdf.pages:
                page_text = page.extract_text()
                if page_text:
                    text_content += page_text + "\n"
            if text_content.strip():
                return text_content
    except Exception as e:
        print(f"pdfplumber failed for {pdf_path}: {str(e)}")

    # Fallback to PyPDF2
    try:
        with open(pdf_path, 'rb') as file:
            pdf_reader = PyPDF2.PdfReader(file)
            text_content = ""
            for page in pdf_reader.pages:
                try:
                    page_text = page.extract_text()
                    if page_text:
                        text_content += page_text + "\n"
                except Exception as e:
                    print(f"Error extracting page from {pdf_path}: {str(e)}")
                    continue
            if text_content.strip():
                return text_content
    except Exception as e:
        print(f"PyPDF2 also failed for {pdf_path}: {str(e)}")

    return ""


def find_best_matching_text(pdf_text, citation, pdf_path):
    """
    Find text in PDF that precisely matches the citation in content and length.

    Args:
        pdf_text (str): Full text content from PDF
        citation (str): Citation text to search for
        pdf_path (Path): Path to PDF file for page-by-page analysis

    Returns:
        tuple: (bool, int, str) - (match_found, page_number, matched_text_from_pdf)
    """
    citation = citation.strip()
    citation_lower = citation.lower()
    citation_length = len(citation)

    # Extract key elements from citation with better precision
    citation_numbers = re.findall(r'\d+\.?\d*', citation)
    citation_units = re.findall(r'\d+\.?\d*\s*[A-Za-z%]+', citation)

    # Get meaningful words, excluding very common ones
    stop_words = {'the', 'and', 'with', 'that', 'were', 'from', 'this', 'has', 'been', 'are', 'was',
                  'for', 'can', 'may', 'also', 'but', 'not', 'more', 'such', 'than', 'into', 'over', 'out'}
    citation_words = [w for w in citation_lower.split() if len(
        w) > 3 and w not in stop_words]

    print(f"    Citation length: {citation_length}")
    print(f"    Looking for numbers: {citation_numbers}")
    print(f"    Looking for units: {citation_units}")
    print(f"    Key words: {citation_words[:7]}")

    # Try pdfplumber first for better text extraction
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page_num, page in enumerate(pdf.pages, 1):
                page_text = page.extract_text()
                if not page_text:
                    continue

                best_match = find_precise_match_in_text(
                    page_text, citation, citation_numbers, citation_units, citation_words, citation_length)
                if best_match:
                    return True, page_num, best_match

    except Exception as e:
        print(f"    pdfplumber failed: {str(e)}")

    # Fallback to PyPDF2
    try:
        with open(pdf_path, 'rb') as file:
            pdf_reader = PyPDF2.PdfReader(file)
            for page_num, page in enumerate(pdf_reader.pages, 1):
                try:
                    page_text = page.extract_text()
                    if not page_text:
                        continue

                    best_match = find_precise_match_in_text(
                        page_text, citation, citation_numbers, citation_units, citation_words, citation_length)
                    if best_match:
                        return True, page_num, best_match

                except Exception as e:
                    continue
    except Exception as e:
        print(f"    PyPDF2 fallback failed: {str(e)}")

    return False, 0, ""


def find_precise_match_in_text(page_text, citation, citation_numbers, citation_units, citation_words, citation_length):
    """
    Find the most precise matching text segment in a page using multiple strategies.
    """
    candidates = []
    citation_lower = citation.lower()
    page_text_lower = page_text.lower()

    # Strategy 1: Direct substring matching (highest priority)
    # Look for exact phrases from citation in PDF
    citation_phrases = []

    # Extract meaningful phrases of different lengths
    words = citation.split()
    for length in [3, 4, 5, 6, 7]:  # Different phrase lengths
        for i in range(len(words) - length + 1):
            phrase = ' '.join(words[i:i+length])
            if len(phrase) > 15:  # Minimum phrase length
                citation_phrases.append(phrase)

    # Check each phrase for exact matches
    for phrase in citation_phrases:
        phrase_lower = phrase.lower()
        if phrase_lower in page_text_lower:
            # Find all occurrences and get context
            start_idx = 0
            while True:
                idx = page_text_lower.find(phrase_lower, start_idx)
                if idx == -1:
                    break

                # Extract context around the match
                context_start = max(0, idx - 100)
                context_end = min(len(page_text), idx + len(phrase) + 100)
                context = page_text[context_start:context_end].strip()

                # Score based on phrase length and position
                score = len(phrase) * 3 + (100 if idx <
                                           len(page_text) // 2 else 50)
                candidates.append(
                    (score, context, f"exact_phrase: {phrase[:30]}..."))

                start_idx = idx + 1

    # Strategy 2: Sentence/line matching with enhanced scoring
    sentences = re.split(r'[.!?;]\s+', page_text)
    lines = [line.strip() for line in page_text.split('\n') if line.strip()]

    for text_segments in [sentences, lines]:
        for segment in text_segments:
            segment = segment.strip()
            if len(segment) < 15:
                continue

            score = calculate_enhanced_score(
                segment, citation, citation_numbers, citation_units, citation_words, citation_length)
            if score > 8:  # Higher threshold
                candidates.append((score, segment, "enhanced_match"))

    # Strategy 3: Multi-line matching for better context (only if no good single matches)
    if not candidates or max(candidates)[0] < 50:
        all_lines = [line.strip()
                     for line in page_text.split('\n') if line.strip()]
        for i in range(len(all_lines)):
            for j in range(i + 1, min(i + 3, len(all_lines) + 1)):  # Max 2 lines
                multi_line = ' '.join(all_lines[i:j])
                if len(multi_line) > citation_length * 1.8:  # Don't get too long
                    break

                score = calculate_enhanced_score(
                    multi_line, citation, citation_numbers, citation_units, citation_words, citation_length)
                if score > 10:
                    candidates.append((score, multi_line, "multi_line"))

    # Return best match
    if candidates:
        candidates.sort(reverse=True, key=lambda x: x[0])
        best_score, best_text, match_type = candidates[0]

        print(f"    Best match score: {best_score}, type: {match_type}")

        # Adaptive threshold based on citation type
        # Lower threshold for numerical citations
        min_threshold = 15 if citation_numbers else 25

        if best_score >= min_threshold:
            return clean_matched_text(best_text)

    return None


def calculate_enhanced_score(text, citation, citation_numbers, citation_units, citation_words, citation_length):
    """
    Enhanced scoring with better phrase and context matching.
    """
    text_lower = text.lower()
    citation_lower = citation.lower()
    score = 0

    # Exact phrase matching (highest priority)
    citation_words_full = citation.split()
    for length in [3, 4, 5]:
        for i in range(len(citation_words_full) - length + 1):
            phrase = ' '.join(citation_words_full[i:i+length]).lower()
            if phrase in text_lower:
                score += len(phrase) * 2

    # Number matching with better context checking
    if citation_numbers:
        text_numbers = re.findall(r'\d+\.?\d*', text)
        for num in citation_numbers:
            if num in text_numbers:
                # Check if number appears with similar units/context
                found_with_context = False
                for unit in citation_units:
                    if num in unit and unit.lower().replace(' ', '') in text_lower.replace(' ', ''):
                        score += 15  # High score for number+unit match
                        found_with_context = True
                        break

                if not found_with_context:
                    score += 5  # Lower score for just number match

    # Key technical terms matching
    technical_keywords = {
        'fracture': 5, 'toughness': 5, 'tensile': 4, 'strength': 4, 'yield': 4,
        'mpa': 6, 'gpa': 6, 'proof': 4, 'stress': 4, 'elongation': 4,
        'oxygen': 5, 'titanium': 5, 'ductility': 5, 'embrittlement': 6,
        'oxide': 4, 'layer': 3, 'powder': 4, 'particles': 3, 'lpbf': 6,
        'alloy': 3, 'properties': 3, 'mechanical': 3
    }

    text_words = text_lower.split()
    for word in citation_words:
        if word in text_words:
            weight = technical_keywords.get(word, 1)
            score += weight

    # Penalty for very different lengths
    length_ratio = abs(len(text) - citation_length) / citation_length
    if length_ratio > 2.0:
        score -= 10
    elif length_ratio < 0.8:
        score += 5

    return score


def clean_matched_text(text):
    """
    Clean and format the matched text.
    """
    # Remove excessive whitespace
    text = re.sub(r'\s+', ' ', text)

    # Remove very short orphaned words at start/end
    words = text.split()
    if len(words) > 3:
        if len(words[0]) <= 2:
            words = words[1:]
        if len(words[-1]) <= 2:
            words = words[:-1]

    return ' '.join(words).strip()


def process_evidence_sources():
    """
    Process evidence sources from assessment.json files and create evidence_source.json
    for each folder containing web search evidence with PDF sources.
    """

    # Base directory path
    SCRIPT_DIR = Path(__file__).parent.absolute()
    base_dir = SCRIPT_DIR / "sonnet45-full-run-codex-gpt5-1010" / "output" / "dry-run"

    # List of folder names to analyze (manually specified)
    folders_to_analyze = [
        "alloys_0003",
        # Add more folder names here as needed
        # "alloys_0014",
        # "alloys_0005",
    ]

    print(
        f"Processing {len(folders_to_analyze)} folders: {folders_to_analyze}\n")

    # Process each folder
    for folder_name in folders_to_analyze:
        folder_path = base_dir / folder_name
        print(f"Processing folder: {folder_name}")

        if not folder_path.exists():
            print(f"Warning: Folder {folder_path} does not exist. Skipping...")
            continue

        # Check for assessment.json
        assessment_path = folder_path / "assessment.json"
        if not assessment_path.exists():
            print(
                f"Warning: assessment.json not found in {folder_path}. Skipping...")
            continue

        # Check for evidences folder
        evidences_path = folder_path / "evidences"
        if not evidences_path.exists():
            print(
                f"Warning: evidences folder not found in {folder_path}. Skipping...")
            continue

        try:
            # Read assessment.json
            with open(assessment_path, 'r', encoding='utf-8') as f:
                assessment_data = json.load(f)

            evidence_mappings = {}

            # Process each evidence item
            if "evidence" in assessment_data:
                for evidence_id, evidence_info in assessment_data["evidence"].items():
                    # Check if it's a web search with URL source
                    if (evidence_info.get("type") == "web search" and
                        "source" in evidence_info and
                            evidence_info["source"].startswith("http")):

                        # Look for corresponding PDF file
                        pdf_filename = f"{evidence_id}.pdf"
                        pdf_path = evidences_path / pdf_filename

                        if pdf_path.exists():
                            print(f"  Processing evidence: {evidence_id}")

                            # Extract PDF text
                            pdf_text = extract_pdf_text(pdf_path)

                            if pdf_text:
                                # Check if citation matches PDF content
                                citation = evidence_info.get("citation", "")
                                if citation:
                                    match_found, page_num, matched_pdf_text = find_best_matching_text(
                                        pdf_text, citation, pdf_path)
                                    if match_found and matched_pdf_text:
                                        # Create evidence mapping in required format
                                        original_url = evidence_info["source"]

                                        evidence_mappings[evidence_id] = {
                                            "original_url": original_url,
                                            "local_path": f"/evidences/{pdf_filename}",
                                            "filename": pdf_filename,
                                            "citation": citation,
                                            "highlights": [
                                                {
                                                    "page": page_num,
                                                    "text": matched_pdf_text,  # Use actual matched PDF text
                                                    "color": "yellow"
                                                }
                                            ]
                                        }
                                        print(
                                            f"    ✓ Added mapping for {evidence_id} (page {page_num})")
                                        print(
                                            f"    Matched text: {matched_pdf_text[:100]}...")
                                    else:
                                        print(
                                            f"    ✗ No matching content found for {evidence_id}")
                                else:
                                    print(
                                        f"    ✗ No citation found for {evidence_id}")
                            else:
                                print(
                                    f"    ✗ Could not extract text from {pdf_filename}")
                        else:
                            print(f"    ✗ PDF file not found: {pdf_filename}")

            # Save evidence_source.json
            if evidence_mappings:
                output_data = {
                    "mappings": evidence_mappings
                }

                output_path = folder_path / "evidence_source.json"
                with open(output_path, 'w', encoding='utf-8') as f:
                    json.dump(output_data, f, indent=4, ensure_ascii=False)

                print(
                    f"  ✓ Saved evidence_source.json with {len(evidence_mappings)} mappings")
            else:
                print(f"  ✗ No evidence mappings found for {folder_name}")

        except Exception as e:
            print(f"Error processing {folder_name}: {str(e)}")
            continue

    print("\nEvidence processing complete!")


if __name__ == "__main__":
    process_evidence_sources()
