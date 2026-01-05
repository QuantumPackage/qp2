#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process IRPF90-generated ninja files to use response files for linking.

This script fixes the "Argument list too long" error by converting linking
commands that have many .o files to use response files (@file syntax).
"""

import re
import sys
import os


def fix_ninja_file(ninja_file_path):
    """
    Read a ninja file and convert long linking commands to use response files.
    
    Args:
        ninja_file_path: Path to the ninja file to fix
    
    Returns:
        True if file was modified, False otherwise
    """
    if not os.path.exists(ninja_file_path):
        return False
    
    with open(ninja_file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    modified = False
    new_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # Look for rule definitions that might be link rules
        if line.startswith('rule '):
            rule_name = line.strip().split()[1] if len(line.strip().split()) > 1 else ''
            rule_lines = [line]
            i += 1
            
            # Collect all lines in this rule
            command_line_text = ''
            while i < len(lines) and (lines[i].startswith('  ') or lines[i].startswith('\t') or lines[i].strip() == ''):
                rule_lines.append(lines[i])
                # Extract command line for keyword checking
                if lines[i].strip().startswith('command ='):
                    command_line_text = lines[i].strip()
                i += 1
            
            # Check if this is a linking rule and needs fixing
            rule_text = ''.join(rule_lines)
            
            # Look for link, archive, or similar commands with $in
            # Check command line specifically for keywords to avoid false positives
            is_link_rule = ('$in' in rule_text and 
                           'rspfile' not in rule_text and
                           command_line_text and
                           any(keyword in command_line_text.lower() 
                               for keyword in ['$fc', '$cc', 'gfortran', 'ifort', 'gcc', 'g++', 'clang', 'ar ']))
            
            if is_link_rule:
                # This rule might benefit from response files
                fixed_rule = fix_link_rule_lines(rule_lines, rule_name)
                new_lines.extend(fixed_rule)
                modified = True
            else:
                new_lines.extend(rule_lines)
        else:
            new_lines.append(line)
            i += 1
    
    if modified:
        # Write the modified content back
        with open(ninja_file_path, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        return True
    
    return False


def fix_link_rule_lines(rule_lines, rule_name):
    """
    Fix a rule to use response files.
    
    Args:
        rule_lines: List of lines making up the rule
        rule_name: Name of the rule
    
    Returns:
        List of fixed lines
    """
    new_rule_lines = []
    command_line_idx = None
    
    # Find the command line
    for idx, line in enumerate(rule_lines):
        if line.strip().startswith('command ='):
            command_line_idx = idx
            break
    
    if command_line_idx is None:
        return rule_lines  # No command found, return as-is
    
    command_line = rule_lines[command_line_idx]
    
    # Replace $in with @$out.rsp in the command using regex
    # Match $in followed by whitespace, end of line, or non-word characters
    new_command = re.sub(r'\$in(?=\s|$|[^\w])', '@$out.rsp', command_line)
    
    # If nothing changed, this rule might not need fixing
    if new_command == command_line:
        return rule_lines
    
    # Build the new rule with rspfile support
    for idx, line in enumerate(rule_lines):
        if idx == command_line_idx:
            new_rule_lines.append(new_command)
            # Add rspfile directives right after the command
            indent = '  ' if line.startswith('  ') else '\t'
            new_rule_lines.append(f'{indent}rspfile = $out.rsp\n')
            new_rule_lines.append(f'{indent}rspfile_content = $in_newline\n')
        else:
            new_rule_lines.append(line)
    
    return new_rule_lines


def process_module_ninja_files(module_path=None):
    """
    Process IRPF90-generated ninja files in a module or the entire source tree.
    
    Args:
        module_path: Path to a specific module directory, or None to process all
    """
    if module_path:
        # Process a specific module
        ninja_file = os.path.join(module_path, 'IRPF90_temp', 'build.ninja')
        if os.path.exists(ninja_file):
            if fix_ninja_file(ninja_file):
                print(f"Fixed: {ninja_file}")
                return 0
            else:
                # File exists but didn't need fixing
                return 0
        else:
            # Ninja file doesn't exist yet, which is fine (might be called too early)
            return 0
    else:
        # Process all modules - for manual use
        qp_root = os.environ.get('QP_ROOT')
        if not qp_root:
            print("Error: QP_ROOT not set", file=sys.stderr)
            return 1
        
        src_dir = os.path.join(qp_root, 'src')
        
        if not os.path.isdir(src_dir):
            print(f"Source directory not found: {src_dir}", file=sys.stderr)
            return 1
        
        count = 0
        for root, dirs, files in os.walk(src_dir):
            if 'IRPF90_temp' in root:
                ninja_file = os.path.join(root, 'build.ninja')
                if os.path.exists(ninja_file):
                    if fix_ninja_file(ninja_file):
                        count += 1
                        print(f"Fixed: {ninja_file}")
        
        if count > 0:
            print(f"\nFixed {count} ninja file(s)")
        else:
            print("No ninja files needed fixing")
        
        return 0


if __name__ == '__main__':
    # Get module path from command line or process all
    if len(sys.argv) > 1:
        module_path = sys.argv[1]
    else:
        module_path = None
    
    sys.exit(process_module_ninja_files(module_path))
