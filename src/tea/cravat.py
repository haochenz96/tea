from pathlib import Path
import time

def write_cravat_input(variants, working_dir, sample_name, *additional_var_info):
    '''
    Writes a cravat input file for a given sample name
    
        inputs:
        - varaints: list of variants to be analyzed (in Tapestri output format)
        - working_dir: directory to place the 'CRAVAT' folder and write cravat input file
        - sample_name: sample name to be used for CRAVAT query

        outputs:

    '''

    # unpack additional_var_info
    additional_var_info = list(additional_var_info)

    if not isinstance(working_dir, Path):
        working_dir = Path(working_dir)

    cravat_input = working_dir / 'CRAVAT' / f'{sample_name}_CRAVAT_input.txt'
    print(f"[INFO] Writing CRAVAT input file to -- {cravat_input}...")
    cravat_input.parents[0].mkdir(
        parents=True, 
        exist_ok=True
    )
    if not cravat_input.is_file():
        with open(cravat_input, 'w') as f:

            var_count=0
            for var in variants:

                var_count += 1
                chrom = var.split(':')[0]
                pos = var.split(':')[1]
                ref = var.split(':')[2].split('/')[0]
                alt = var.split(':')[2].split('/')[1]

                f.write(f'{chrom} {pos} + {ref} {alt} var_{var_count} {" ".join([str(additional_var_info[i][var]) for i in range(len(additional_var_info))])}\n')
        print(f"[INFO] Finished writing CRAVAT input file!")
    else:
        print(f"[WARNING] CRAVAT input file already exists!")

def get_cravat_output(session, job_id, output_path):
    '''
    Writes a cravat input file for a given sample name
    
        inputs:
        - session: requests.session where the cravat query is in
        - job_id: can be feteched by post.json()['id']
        - output_path: path to write the CRAVAT output file

        outputs:

    '''

    response = session.get(
        'https://run.opencravat.org/submit/jobs/' + job_id + '/status')
    if response.json()['status'] == 'error':
        print(f'[WARNING] CRAVAT run failed! Check input file')
        
    elif response.json()['status'] != 'Finished':
        time.sleep(5)
        get_cravat_output(session, job_id, output_path)

    elif response.json()['status'] == 'Finished':
        output = session.get('https://run.opencravat.org/submit/jobs/' + job_id + '/reports/text')

        with open(output_path, 'w') as f:
            f.write(output.text)

    else:
        print(f'[WARNING] CRAVAT run failed! Unknown error')