#!/usr/bin/env python3
"""
Script to fetch biosample IDs associated with a bioproject ID from NCBI.

Requirements:
    pip install requests

Usage:
    python getBiosamplefromProjectID.py --bioproject PRJNA123456 --output results.txt
    python getBiosamplefromProjectID.py -b PRJNA123456 -o results.csv --format csv
    python getBiosamplefromProjectID.py PRJNA123456  # Simple usage (saves as PRJNA123456_biosamples.txt)
"""

import requests
import time
import argparse
import sys
from typing import List, Optional

class NCBIBioprojectFetcher:
    def __init__(self, email: Optional[str] = None):
        """
        Initialize the NCBI fetcher.
        
        Args:
            email: Your email address (recommended by NCBI for API usage)
        """
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.email = email
        
    def fetch_biosamples(self, bioproject_id: str) -> List[str]:
        """
        Fetch biosample IDs from NCBI for a given bioproject ID.
        
        Args:
            bioproject_id: NCBI bioproject ID (e.g., 'PRJNA123456')
            
        Returns:
            List of biosample accession IDs
        """
        print(f"Searching NCBI for bioproject: {bioproject_id}")
        
        try:
            # Step 1: Search for the bioproject to get its UID
            bioproject_uid = self._search_bioproject(bioproject_id)
            if not bioproject_uid:
                return []
            
            # Step 2: Find linked biosample UIDs
            biosample_uids = self._get_linked_biosamples(bioproject_uid)
            if not biosample_uids:
                return []
            
            # Step 3: Convert UIDs to accession IDs
            biosample_ids = self._get_biosample_accessions(biosample_uids)
            
            return biosample_ids
            
        except requests.exceptions.RequestException as e:
            print(f"Network error: {e}")
            return []
        except Exception as e:
            print(f"Unexpected error: {e}")
            return []
    
    def _search_bioproject(self, bioproject_id: str) -> Optional[str]:
        """Search for bioproject and return its UID."""
        params = {
            'db': 'bioproject',
            'term': bioproject_id,
            'retmode': 'json'
        }
        
        if self.email:
            params['email'] = self.email
        
        response = requests.get(f"{self.base_url}esearch.fcgi", params=params)
        response.raise_for_status()
        
        data = response.json()
        
        if not data['esearchresult']['idlist']:
            print(f"No bioproject found for ID: {bioproject_id}")
            return None
        
        uid = data['esearchresult']['idlist'][0]
        print(f"Found bioproject UID: {uid}")
        return uid
    
    def _get_linked_biosamples(self, bioproject_uid: str) -> List[str]:
        """Get biosample UIDs linked to the bioproject."""
        params = {
            'dbfrom': 'bioproject',
            'db': 'biosample',
            'id': bioproject_uid,
            'retmode': 'json'
        }
        
        if self.email:
            params['email'] = self.email
        
        response = requests.get(f"{self.base_url}elink.fcgi", params=params)
        response.raise_for_status()
        
        data = response.json()
        
        biosample_uids = []
        if 'linksets' in data and data['linksets']:
            for linkset in data['linksets']:
                if 'linksetdbs' in linkset:
                    for linksetdb in linkset['linksetdbs']:
                        if linksetdb['dbto'] == 'biosample':
                            biosample_uids.extend(linksetdb['links'])
        
        if not biosample_uids:
            print("No linked biosamples found")
            return []
        
        print(f"Found {len(biosample_uids)} linked biosample UIDs")
        return biosample_uids
    
    def _get_biosample_accessions(self, biosample_uids: List[str]) -> List[str]:
        """Convert biosample UIDs to accession IDs."""
        # Process UIDs in batches to avoid URL length limits
        batch_size = 200
        all_accessions = []
        
        for i in range(0, len(biosample_uids), batch_size):
            batch_uids = biosample_uids[i:i + batch_size]
            
            params = {
                'db': 'biosample',
                'id': ','.join(batch_uids),
                'retmode': 'json'
            }
            
            if self.email:
                params['email'] = self.email
            
            response = requests.get(f"{self.base_url}esummary.fcgi", params=params)
            response.raise_for_status()
            
            data = response.json()
            
            if 'result' in data:
                for uid in batch_uids:
                    if uid in data['result']:
                        accession = data['result'][uid].get('accession', '')
                        if accession:
                            all_accessions.append(accession)
            
            # Be respectful to NCBI servers
            if i + batch_size < len(biosample_uids):
                time.sleep(0.5)
        
        return sorted(all_accessions)
    
    def save_to_file(self, biosamples: List[str], bioproject_id: str, filename: Optional[str] = None) -> str:
        """
        Save biosample IDs to a file.
        
        Args:
            biosamples: List of biosample IDs
            bioproject_id: The bioproject ID
            filename: Optional custom filename
            
        Returns:
            The filename where data was saved
        """
        if filename is None:
            filename = f"{bioproject_id}_biosamples.txt"
        
        try:
            with open(filename, 'w') as f:
                f.write(f"Biosample IDs for bioproject: {bioproject_id}\n")
                f.write(f"Retrieved from NCBI on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total count: {len(biosamples)}\n")
                f.write("=" * 50 + "\n\n")
                
                for i, biosample_id in enumerate(biosamples, 1):
                    f.write(f"{i:4d}. {biosample_id}\n")
                
                f.write(f"\n" + "=" * 50 + "\n")
                f.write("End of file\n")
            
            print(f"Results saved to: {filename}")
            return filename
            
        except IOError as e:
            print(f"Error saving file: {e}")
            return ""
    
    def save_to_csv(self, biosamples: List[str], bioproject_id: str, filename: Optional[str] = None) -> str:
        """
        Save biosample IDs to a CSV file.
        
        Args:
            biosamples: List of biosample IDs
            bioproject_id: The bioproject ID
            filename: Optional custom filename
            
        Returns:
            The filename where data was saved
        """
        if filename is None:
            filename = f"{bioproject_id}_biosamples.csv"
        
        try:
            with open(filename, 'w') as f:
                f.write("bioproject_id,biosample_id,index\n")
                for i, biosample_id in enumerate(biosamples, 1):
                    f.write(f"{bioproject_id},{biosample_id},{i}\n")
            
            print(f"CSV results saved to: {filename}")
            return filename
            
        except IOError as e:
            print(f"Error saving CSV file: {e}")
            return ""
    
    def fetch_and_save(self, bioproject_id: str, save_format: str = 'txt', filename: Optional[str] = None) -> dict:
        """
        Fetch biosamples and automatically save to file.
        
        Args:
            bioproject_id: NCBI bioproject ID
            save_format: 'txt' or 'csv'
            filename: Optional custom filename
            
        Returns:
            Dictionary with results and saved filename
        """
        biosamples = self.fetch_biosamples(bioproject_id)
        
        result = {
            'bioproject_id': bioproject_id,
            'biosample_count': len(biosamples),
            'biosamples': biosamples,
            'saved_file': ''
        }
        
        if biosamples:
            if save_format.lower() == 'csv':
                result['saved_file'] = self.save_to_csv(biosamples, bioproject_id, filename)
            else:
                result['saved_file'] = self.save_to_file(biosamples, bioproject_id, filename)
        
        return result
    
    def fetch_with_details(self, bioproject_id: str) -> dict:
        """
        Fetch biosamples with additional details.
        
        Args:
            bioproject_id: NCBI bioproject ID
            
        Returns:
            Dictionary with biosample details
        """
        biosample_ids = self.fetch_biosamples(bioproject_id)
        
        if not biosample_ids:
            return {'bioproject_id': bioproject_id, 'biosample_count': 0, 'biosamples': []}
        
        return {
            'bioproject_id': bioproject_id,
            'biosample_count': len(biosample_ids),
            'biosamples': biosample_ids
        }

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Fetch biosample IDs from NCBI for a given bioproject ID',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --bioproject PRJNA123456 --output results.txt
  %(prog)s -b PRJNA123456 -o results.csv --format csv
  %(prog)s PRJNA123456  # Simple usage, saves as PRJNA123456_biosamples.txt
  %(prog)s -b PRJNA123456 --email user@example.com --format csv
        """
    )
    
    # Positional argument for bioproject ID (optional if --bioproject is used)
    parser.add_argument('bioproject_id', nargs='?', help='Bioproject ID (e.g., PRJNA123456)')
    
    # Optional arguments
    parser.add_argument('-b', '--bioproject', help='Bioproject ID (alternative to positional argument)')
    parser.add_argument('-o', '--output', help='Output filename (default: auto-generated based on bioproject ID)')
    parser.add_argument('-f', '--format', choices=['txt', 'csv'], default='txt', help='Output format (default: txt)')
    parser.add_argument('-e', '--email', help='Email address (recommended by NCBI)')
    parser.add_argument('--no-save', action='store_true', help='Print results only, do not save to file')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress progress messages')
    parser.add_argument('--list-only', action='store_true', help='Output only the biosample IDs, one per line')
    
    args = parser.parse_args()
    
    # Determine bioproject ID from positional or named argument
    bioproject_id = args.bioproject_id or args.bioproject
    if not bioproject_id:
        parser.error("Bioproject ID is required. Use positional argument or --bioproject/-b option.")
    
    args.bioproject_id = bioproject_id
    return args

def main():
    """Main function with command line argument support."""
    args = parse_arguments()
    
    # Initialize fetcher
    fetcher = NCBIBioprojectFetcher(email=args.email)
    
    if not args.quiet:
        print(f"Fetching biosample IDs for bioproject: {args.bioproject_id}")
        print("-" * 50)
    
    # Fetch biosamples
    biosamples = fetcher.fetch_biosamples(args.bioproject_id)
    
    # Handle results
    if not biosamples:
        if not args.quiet:
            print("No biosample IDs found.")
        sys.exit(1)
    
    if args.list_only:
        # Just print the list, one per line
        for biosample in biosamples:
            print(biosample)
        return
    
    if not args.quiet:
        print(f"Found {len(biosamples)} biosample IDs:")
        print("-" * 30)
        for i, biosample_id in enumerate(biosamples, 1):
            print(f"{i:4d}. {biosample_id}")
    
    # Save results unless --no-save is specified
    if not args.no_save:
        # Determine output filename
        if args.output:
            output_file = args.output
        else:
            extension = 'csv' if args.format == 'csv' else 'txt'
            output_file = f"{args.bioproject_id}_biosamples.{extension}"
        
        # Save based on format
        if args.format == 'csv':
            saved_file = fetcher.save_to_csv(biosamples, args.bioproject_id, output_file)
        else:
            saved_file = fetcher.save_to_file(biosamples, args.bioproject_id, output_file)
        
        if not args.quiet and saved_file:
            print(f"\nResults saved to: {saved_file}")

def interactive_main():
    """Original interactive main function for backward compatibility."""
    print("NCBI Bioproject to Biosample Fetcher")
    print("=" * 40)
    
    # Get user input
    bioproject_id = input("Enter Bioproject ID (e.g., PRJNA123456): ").strip()
    email = input("Enter your email (optional but recommended): ").strip()
    
    if not bioproject_id:
        print("No bioproject ID provided. Using example: PRJNA123456")
        bioproject_id = "PRJNA123456"
    
    # Initialize fetcher
    fetcher = NCBIBioprojectFetcher(email=email if email else None)
    
    print(f"\nFetching biosample IDs for bioproject: {bioproject_id}")
    print("-" * 50)
    
    # Fetch biosamples
    biosamples = fetcher.fetch_biosamples(bioproject_id)
    
    # Display results
    print(f"\nResults:")
    print(f"Found {len(biosamples)} biosample IDs:")
    print("-" * 30)
    
    if biosamples:
        for i, biosample_id in enumerate(biosamples, 1):
            print(f"{i:4d}. {biosample_id}")
        
        # Save to file options
        print(f"\nSave options:")
        print("1. Save as text file")
        print("2. Save as CSV file") 
        print("3. Don't save")
        
        save_choice = input("Choose option (1-3): ").strip()
        
        if save_choice == "1":
            filename = f"{bioproject_id}_biosamples.txt"
            saved_file = fetcher.save_to_file(biosamples, bioproject_id)
        elif save_choice == "2":
            filename = f"{bioproject_id}_biosamples.csv"
            saved_file = fetcher.save_to_csv(biosamples, bioproject_id)
        else:
            print("Results not saved.")
    else:
        print("No biosample IDs found.")

# Example usage as a module
def example_usage():
    """Example of how to use this as a module."""
    # Initialize fetcher
    fetcher = NCBIBioprojectFetcher(email="your.email@example.com")
    
    # Method 1: Fetch and save automatically
    print("Method 1: Fetch and save automatically")
    result = fetcher.fetch_and_save("PRJNA123456", save_format='txt')
    print(f"Found {result['biosample_count']} biosamples, saved to: {result['saved_file']}")
    
    # Method 2: Fetch then save manually
    print("\nMethod 2: Fetch then save manually")
    biosamples = fetcher.fetch_biosamples("PRJNA123456")
    if biosamples:
        # Save as text file
        fetcher.save_to_file(biosamples, "PRJNA123456", "my_custom_filename.txt")
        
        # Save as CSV file
        fetcher.save_to_csv(biosamples, "PRJNA123456", "my_custom_filename.csv")
    
    # Method 3: Just get the list (no saving)
    print("\nMethod 3: Just get the list")
    biosamples = fetcher.fetch_biosamples("PRJNA123456")
    print(f"Found {len(biosamples)} biosamples: {biosamples[:5]}...")  # Show first 5

if __name__ == "__main__":
    # Check if any command line arguments were provided
    if len(sys.argv) > 1:
        main()  # Use command line interface
    else:
        interactive_main()  # Fall back to interactive mode